"""Mash sequence neighborhoods and near-clone components, with streaming edges."""

import gzip
import json
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .audit import safe_path, sha256, write_tsv
from .cohort import UnionFind, read_table
from .snapshot import atomic_json, prepare_output


def run_diversity(representatives_path, config_path, output, threads=8):
    representatives_path, config_path, output = map(Path, (representatives_path, config_path, output))
    config = json.loads(config_path.read_text())
    if not 0 <= config["near_clone_distance"] <= config["neighborhood_distance"] <= 1:
        raise ValueError("Require 0 <= near_clone_distance <= neighborhood_distance <= 1")
    if config["sketch_batch_size"] < 1:
        raise ValueError("sketch_batch_size must be positive")
    if "raw" in output.resolve().parts or not 1 <= threads <= 24:
        raise ValueError("Invalid output or threads")
    records = read_table(representatives_path)
    if not records:
        raise ValueError("No eligible representatives")
    binary = Path(config["mash_binary"])
    version = subprocess.check_output([str(binary), "--version"], text=True).strip()
    root = Path(config["genome_root"])
    names = {str(safe_path(root, r["target_relpath"])): r["genome_id"] for r in records}
    if len(names) != len(records):
        raise ValueError("Duplicate sequence paths")
    contract = {"representatives_sha256": sha256(representatives_path), "config_sha256": sha256(config_path), "mash_binary_sha256": sha256(binary), "mash_version": version, "code_sha256": sha256(__file__), "threads": threads}
    prepare_output(output, contract)
    def verify_inputs():
        for r in records:
            path = safe_path(root, r["target_relpath"])
            stat = path.stat()
            if (stat.st_size, stat.st_mtime_ns) != (int(r["initial_size"]), int(r["initial_mtime_ns"])):
                raise ValueError(f"Input changed since snapshot: {r['genome_id']}")
            if sha256(path) != r["file_sha256"]:
                raise ValueError(f"Input content changed since snapshot: {r['genome_id']}")
    verify_inputs()
    completion = output / "completion.json"
    if completion.exists():
        for name, digest in json.loads(completion.read_text()).items():
            if sha256(safe_path(output, name)) != digest:
                raise ValueError(f"Completed artifact changed: {name}")
        atomic_json(output / "progress.json", {"status": "completed", "records": len(records)})
        return json.loads((output / "summary.json").read_text())
    sketches = []
    paths = list(names)
    for start in range(0, len(paths), config["sketch_batch_size"]):
        prefix = output / f"sketch_{start:06d}"
        sketch = prefix.with_suffix(".msh")
        done = prefix.with_suffix(".done.json")
        inputs = prefix.with_suffix(".inputs.txt")
        inputs.write_text("\n".join(paths[start:start + config["sketch_batch_size"]]) + "\n")
        if not (sketch.exists() and done.exists() and json.loads(done.read_text())["sha256"] == sha256(sketch)):
            command = [str(binary), "sketch", "-l", "-p", str(threads), "-k", str(config["k"]), "-s", str(config["sketch_size"]), "-S", str(config["seed"]), "-o", str(prefix), str(inputs)]
            with prefix.with_suffix(".log").open("w") as log:
                subprocess.run(command, check=True, stdout=log, stderr=log)
            atomic_json(done, {"command": command, "sha256": sha256(sketch)})
        sketches.append(str(sketch))
        atomic_json(output / "progress.json", {"status": "sketching", "completed": min(start + config["sketch_batch_size"], len(paths)), "total": len(paths)})
        print(f"sketch {min(start + config['sketch_batch_size'], len(paths))}/{len(paths)}", flush=True)
    combined = output / "genomes.msh"
    temporary = output / "genomes.pending.msh"
    temporary.unlink(missing_ok=True)
    subprocess.run([str(binary), "paste", str(output / "genomes.pending"), *sketches], check=True, stdout=subprocess.DEVNULL)
    temporary.replace(combined)
    info = subprocess.check_output([str(binary), "info", "-t", str(combined)], text=True)
    (output / "sketch_inventory.tsv").write_text(info)
    observed = [line.split("\t")[2] for line in info.splitlines() if line and not line.startswith("#Hashes\t")]
    if len(observed) != len(records) or set(observed) != set(names):
        raise ValueError("Mash inventory does not match cohort")
    near = UnionFind(names.values())
    neighborhoods = UnionFind(names.values())
    command = [str(binary), "triangle", "-E", "-d", str(config["neighborhood_distance"]), "-p", str(threads), str(combined)]
    atomic_json(output / "progress.json", {"status": "pairwise_distances", "total_pairs": len(records) * (len(records) - 1) // 2})
    edge_count, near_count, forest_count, broad_forest_count = 0, 0, 0, 0
    with (output / "triangle.log").open("w") as log, gzip.open(output / "near_clone_forest.tsv.gz", "wt") as edges, gzip.open(output / "neighborhood_forest.tsv.gz", "wt") as broad_edges:
        edges.write("genome_a\tgenome_b\tmash_distance\tshared_hashes\n")
        broad_edges.write("genome_a\tgenome_b\tmash_distance\tshared_hashes\n")
        with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log, text=True, bufsize=1) as process:
            for line in process.stdout:
                a, b, distance, pvalue, shared = line.rstrip().split("\t")
                a, b, distance = names[a], names[b], float(distance)
                if distance > config["neighborhood_distance"]:
                    raise ValueError("Mash returned an edge outside requested distance")
                if neighborhoods.union(a, b):
                    broad_forest_count += 1
                    broad_edges.write(f"{a}\t{b}\t{distance}\t{shared}\n")
                edge_count += 1
                if distance <= config["near_clone_distance"]:
                    near_count += 1
                    # A spanning forest proves exactly the same components while
                    # avoiding quadratic storage for dense clone neighborhoods.
                    if near.union(a, b):
                        forest_count += 1
                        edges.write(f"{a}\t{b}\t{distance}\t{shared}\n")
                if edge_count % 1000000 == 0:
                    atomic_json(output / "progress.json", {"status": "pairwise_distances", "edges_within_neighborhood_threshold": edge_count, "near_clone_edges": near_count})
                    print(f"distance edges {edge_count}, near-clone edges {near_count}", flush=True)
            if process.wait() != 0:
                raise RuntimeError("Mash triangle failed; see triangle.log")
    verify_inputs()
    assignments = [{"genome_id": gid, "near_clone_group": near.find(gid), "sequence_neighborhood": neighborhoods.find(gid)} for gid in sorted(names.values())]
    write_tsv(output / "sequence_groups.tsv", assignments, ["genome_id", "near_clone_group", "sequence_neighborhood"])
    near_groups, broad_groups = near.groups(), neighborhoods.groups()
    if forest_count != len(records) - len(near_groups) or broad_forest_count != len(records) - len(broad_groups):
        raise ValueError("Stored forests do not cover the sequence components")
    summary = {"status": "sequence_grouping_complete", "records": len(records), "near_clone_groups": len(near_groups), "sequence_neighborhoods": len(broad_groups), "near_clone_component_sizes": dict(Counter(len(g) for g in near_groups.values())), "neighborhood_component_sizes": dict(Counter(len(g) for g in broad_groups.values())), "near_clone_edges": near_count, "neighborhood_edges": edge_count, "config": config, "command": command, "finished_utc": datetime.now(timezone.utc).isoformat(), "groups_sha256": sha256(output / "sequence_groups.tsv"), "sketch_sha256": sha256(combined), "limitations": ["Mash distances are MinHash estimates of whole-assembly k-mer divergence, not core SNP counts or exact ANI.", "Near-clone and neighborhood groups are single-linkage connected components; not all pairs within a component meet the edge threshold.", "Sequence neighborhoods are operational genetic structure, not validated MLST or phylogenetic lineages.", "Sketch resolution and accessory content affect grouping; source/BioSample identity is handled separately."]}
    summary.update(near_clone_forest_edges=forest_count, neighborhood_forest_edges=broad_forest_count)
    atomic_json(output / "summary.json", summary)
    atomic_json(completion, {name: sha256(output / name) for name in ("summary.json", "sequence_groups.tsv", "genomes.msh", "sketch_inventory.tsv", "near_clone_forest.tsv.gz", "neighborhood_forest.tsv.gz")})
    atomic_json(output / "progress.json", {"status": "completed", "records": len(records)})
    return json.loads((output / "summary.json").read_text())
