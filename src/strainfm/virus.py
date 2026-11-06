"""Split a virus archive into a new directory only after complete decompression."""

import csv
import gzip
import re
import subprocess
import tempfile
from pathlib import Path

from .audit import safe_path


def split_virus(manifest, archive, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError(f"Refusing to replace existing raw output: {output}")
    with Path(manifest).open(newline="", encoding="utf-8") as handle:
        mapping = {r["genome_id"]: r["species_id"] for r in csv.DictReader(handle, delimiter="\t") if r["genome_id"].startswith("UHGV-")}
    output.parent.mkdir(parents=True, exist_ok=True)
    written, unparsed, unmapped = set(), 0, 0
    with tempfile.TemporaryDirectory(dir=output.parent, prefix=".uhgv-") as directory:
        stage = Path(directory) / "genomes"
        stage.mkdir()
        current = None
        with subprocess.Popen(["zstd", "--long=29", "--memory=1024MB", "-dc", str(archive)], stdout=subprocess.PIPE, text=True, encoding="utf-8") as process:
            try:
                for line in process.stdout:
                    if line.startswith(">"):
                        if current is not None:
                            current.close()
                            current = None
                        match = re.search(r"UHGV-\d+", line)
                        if not match:
                            unparsed += 1
                            continue
                        gid = match.group(0)
                        if gid not in mapping:
                            unmapped += 1
                            continue
                        if gid in written:
                            raise ValueError(f"Duplicate virus record: {gid}")
                        target = safe_path(stage, str(Path(mapping[gid]) / (gid + ".fna.gz")))
                        target.parent.mkdir(parents=True, exist_ok=True)
                        current = gzip.open(target, "xt", encoding="utf-8")
                        current.write(line)
                        written.add(gid)
                    elif current is not None:
                        current.write(line)
                if process.wait() != 0:
                    raise RuntimeError("Virus decompression failed; no output published")
            except BaseException:
                if process.poll() is None:
                    process.kill()
                process.wait()
                raise
            finally:
                if current is not None:
                    current.close()
        if output.exists():
            raise FileExistsError(output)
        stage.rename(output)
    return {"records": len(written), "manifest_missing": len(set(mapping) - written), "unparsed_headers": unparsed, "unmapped_headers": unmapped}
