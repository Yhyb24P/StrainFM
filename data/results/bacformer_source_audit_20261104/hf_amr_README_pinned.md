---
license: apache-2.0
tags:
- antibiotic
- bacteria
- genomics
- proteins
- biology
- prokaryotes
pretty_name: >-
  Dataset for predicting antibiotic resistance from whole bacterial genomes
  (protein sequences)
size_categories:
- 10K<n<100K
---
# Dataset for antibiotic resistance prediction from whole-bacterial genomes (protein sequences)

A dataset of 25,032 bacterial genomes across 39 species with antimicrobial resistance labels.

The genome protein sequences have been extracted from [GenBank](https://www.ncbi.nlm.nih.gov/genbank/). Each row contains whole bacterial genome, with spaces
separating different contigs present in the genome.

The antimicrobial resistance labels have been extracted from [Antibiotic Susceptibility Test (AST) Browser](https://www.ncbi.nlm.nih.gov/pathogens/ast), accessed 23 Oct, 2024.)
and include both `binary` (resistant/susceptible) labels as well as `minimum inhibitory concentration (MIC)` regression values. The MIC has been `log1p` normalised.

We exclude antimicrobials with a low nr of samples, giving us `56` unique antimicrobials for `MIC (regression)` prediction and `36` for binary labels. For binary case, we only
included genomes with `susceptible` and `resistant` labels provided by the AST Browser, excluding ambiguous labels. We treat combination of antimicrobials as a separate drug.

## Labels
We provide labels in separate files in the dataset `Files and versions`. This includes:
* binary labels - `binary_labels.csv`
* MIC (regression) labels - `mic_regression_labels.csv`

## Usage
We recommend loading the dataset in a streaming mode to prevent memory errors.
```python
from datasets import load_dataset


ds = load_dataset("macwiatrak/bacbench-antibiotic-resistance-dna", split="train", streaming=True)
```

### Fetch the labels for the genome
```python
import pandas as pd

from datasets import load_dataset


ds = load_dataset("macwiatrak/bacbench-antibiotic-resistance-protein-sequences", split="train", streaming=True)
item = next(iter(ds))

# read labels (available in repo root)
labels_df = pd.read_csv("<input-dir>/mic_regression_labels.csv").set_index("genome_name")

# fetch labels
labels = labels_df.loc[item["genome_name"]]
# drop antibiotics without a value for the genome (NaN)
labels = labels.dropna()
```

## Split

Due to low number of samples for many antibiotics and the variability between genomes, which may skew the results when using a single split, we recommend training and evaluating the model with `k-fold split`. 
Specifically, for each antibiotic we recommend:
1. Splitting the available data into 5 equal splits (`sklearn.model_selection.StratifiedKFold` for binary labels and `sklearn.model_selection.KFold` for regression labels)
2. In each split, further dividing the larger `train` set into `train` and `val`, where `validation` makes up 20% of the train split.
3. Training the model on the train set from the point above and monitoring the results on the validation set, using `AUPRC` and `R2` as metrics for monitoring the performance on the validation set for binary and regression setups.
4. Using the best performing model on validation to evaluate the model on the test set.

**NOTE**: If your model is going to be trained on multiple antibiotics concurrently, the same genome will appear in the set multiple times. It is then important to split to create train/validation/test sets on genome_id (and expand the splits to include all tests on those genomes) to prevent genome leak into the test set.

See [github repository](https://github.com/macwiatrak/Bacbench) for details on how to embed the dataset with DNA and protein language models as well as code to predict antibiotic
resistance from sequence. For DNA sequence representation of the genomes see the [antibiotic-resistance-protein-sequences](https://huggingface.co/datasets/macwiatrak/bacbench-antibiotic-resistance-dna)
dataset.

---
dataset_info:
  features:
  - name: genome_name
    dtype: string
  - name: contig_name
    sequence: string
  - name: protein_id
    sequence:
      sequence: string
  - name: protein_sequence
    sequence:
      sequence: string
  - name: taxid
    dtype: string
  - name: locus_tag
    sequence:
      sequence: string
  - name: start
    sequence:
      sequence: int64
  - name: end
    sequence:
      sequence: int64
  - name: product
    sequence:
      sequence: string
  splits:
  - name: train
    num_bytes: 40490436864
    num_examples: 26052
  download_size: 34207458365
  dataset_size: 40490436864
configs:
- config_name: default
  data_files:
  - split: train
    path: data/train-*
license: apache-2.0
tags:
- biology
- antibiotic
- resistance
- bacteria
- genomics
size_categories:
- 10K<n<100K
---