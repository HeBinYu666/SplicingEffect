# IsoImpact

IsoImpact is an automated pipeline designed to evaluate the functional impact of alternative splicing. It systematically compares protein sequence features, genomic coordinates, and domain annotations across different transcript isoforms to reveal the functional consequences of alternative splicing events.

## 1. Installation

Download IsoImpact from GitHub:

```bash
git clone https://github.com/HeBinYu666/IsoImpact.git
cd IsoImpact
```

IsoImpact is implemented in Python. It requires Python 3.8 or higher. The required Python packages are listed in `requirements.txt`:

```text
numpy>=1.21
pandas>=1.3
matplotlib>=3.5
biopython>=1.79
propy3>=1.1.1
```

Install these packages with:

```bash
pip install -r requirements.txt
```

No R package is required for the standard human and mouse IsoImpact workflows. R is only needed when users want to build a custom domain-coordinate file for novel or custom isoforms.A dedicated tutorial for this process is provided in Part 3.

## 2. Usage

IsoImpact requires three main inputs:

```text
1. Ensembl transcript IDs for the isoforms to be compared.
2. A matching GTF annotation file.
3. A matching protein FASTA file.
```

The transcript IDs must match the supplied GTF file. The GTF file is required because IsoImpact uses it to identify transcript, exon, CDS, UTR, gene, and protein ID annotations for the selected isoforms. The protein FASTA file is required because IsoImpact extracts the corresponding amino acid sequences and calculates protein sequence features from those sequences.

For standard human and mouse workflows, the GTF and protein FASTA files should come from Ensembl release 110, because the built-in domain-coordinate files in this repository were generated from the same Ensembl release. Using matching files helps ensure that transcript IDs, protein IDs, genomic coordinates, and domain annotations can be linked correctly. These reference GTF and protein FASTA files are large, so they are not included directly in this GitHub repository. Users should download them from Ensembl before running IsoImpact.

The main command uses the following options:

```text
-i / --isoform   Transcript IDs to compare. At least two transcript IDs are required.
-g / --gtf       Path to the matching GTF annotation file.
-f / --fasta     Path to the matching protein FASTA file.
-o / --outdir    Directory where output files will be saved. The default is the current directory.
-d / --domain    Optional custom domain-coordinate CSV file. This is not needed for standard human or mouse Ensembl release 110 runs.
```

For example, if `-o results/human_isoforms` is used, IsoImpact writes the main output files as:

```text
results/human_isoforms/IsoImpact_features.csv
results/human_isoforms/IsoImpact_figure.png
```

### 2.1 Human

Download the Ensembl release 110 human GTF and protein FASTA files:

```bash
mkdir -p references

curl -L -o references/Homo_sapiens.GRCh38.110.gtf.gz \
  https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

curl -L -o references/Homo_sapiens.GRCh38.pep.all.fa.gz \
  https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz

gunzip references/Homo_sapiens.GRCh38.110.gtf.gz
gunzip references/Homo_sapiens.GRCh38.pep.all.fa.gz
```

Run IsoImpact with human Ensembl transcript IDs:

```bash
python IsoImpact.py \
  -i ENST00000421030 ENST00000440047 ENST00000413188 ENST00000409996 \
  -g references/Homo_sapiens.GRCh38.110.gtf \
  -f references/Homo_sapiens.GRCh38.pep.all.fa \
  -o results/human_isoforms
```

For human transcript IDs beginning with `ENST`, IsoImpact automatically uses:

```text
data/human_domain.csv
```

Therefore, `-d/--domain` is not required for the standard human Ensembl release 110 workflow.

### 2.2 Mouse

Download the Ensembl release 110 mouse GTF and protein FASTA files:

```bash
mkdir -p references

curl -L -o references/Mus_musculus.GRCm39.110.gtf.gz \
  https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz

curl -L -o references/Mus_musculus.GRCm39.pep.all.fa.gz \
  https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz

gunzip references/Mus_musculus.GRCm39.110.gtf.gz
gunzip references/Mus_musculus.GRCm39.pep.all.fa.gz
```

Run IsoImpact with mouse Ensembl transcript IDs:

```bash
python IsoImpact.py \
  -i ENSMUST00000193361 ENSMUST00000064562 ENSMUST00000092420 ENSMUST00000105424 \
  -g references/Mus_musculus.GRCm39.110.gtf \
  -f references/Mus_musculus.GRCm39.pep.all.fa \
  -o results/mouse_isoforms
```

Replace the example transcript IDs with the mouse isoforms to be compared. For mouse transcript IDs beginning with `ENSMUST`, IsoImpact automatically uses:

```text
data/mouse_domain.csv
```

Therefore, `-d/--domain` is not required for the standard mouse Ensembl release 110 workflow.

## 3. Custom Domain Annotation

This section is only needed for novel isoforms, custom isoforms, non-model organisms, or annotations that are not covered by the built-in human and mouse Ensembl release 110 domain-coordinate files.

For a custom analysis, users need:

```text
1. A CDS-aware GTF file.
2. A matching protein FASTA file.
3. Pfam/PfamScan domain-prediction results for the corresponding protein sequences(e.g., generated via local PfamScan or InterProScan).
```

The custom GTF file must contain CDS records, because the helper script maps protein-domain intervals back to genomic CDS coordinates. The protein IDs in the Pfam/PfamScan result file must match the `protein_id` values in the GTF file. The same protein IDs should also be used in the protein FASTA headers.

Install the R packages required by the helper script:

```r
install.packages("BiocManager")
BiocManager::install(c("ensembldb", "GenomicRanges"))
install.packages("dplyr")
```

Build an IsoImpact-compatible custom domain-coordinate file:

```bash
Rscript scripts/build_custom_domain.R \
  --gtf your_novel_isoforms.gtf \
  --pfam your_pfam_results.txt \
  --out custom_domain.csv
```

Then run IsoImpact with the custom domain-coordinate file:

```bash
python IsoImpact.py \
  -i novel_tx_1 novel_tx_2 \
  -g your_novel_isoforms.gtf \
  -f your_novel_proteins.fa \
  -d custom_domain.csv \
  -o results/novel
```

Use `-d/--domain` only when providing a custom domain-coordinate file. Standard human and mouse Ensembl release 110 workflows do not require this option.

## 4. Output Files

IsoImpact generates two main output files:

```text
IsoImpact_features.csv
IsoImpact_figure.png
```

`IsoImpact_features.csv` contains:

```text
gene, transcript, protein, and coding biotype metadata
protein length and molecular weight comparisons
exon count and genomic span comparisons
shared, lost, and gained domain summaries
canonical, alternative, and delta values for protein sequence features
```

`IsoImpact_figure.png` contains the domain-mapping visualization and ranked feature-difference plots for the input isoforms.

IsoImpact selects the canonical baseline by parsing `Ensembl_canonical` or `MANE_Select` tags when available. If neither tag is available, it uses the input transcript with the largest genomic span.

## 5. Contact

Please open a GitHub issue for questions, bug reports, or feature requests.
