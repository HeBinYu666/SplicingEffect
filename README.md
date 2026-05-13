# IsoImpact

## 1. Installation

Clone the repository and enter the project directory:

```bash
git clone https://github.com/HeBinYu666/IsoImpact.git
cd IsoImpact
```

Create and activate a Python environment:

```bash
python -m venv .venv
source .venv/bin/activate
```

On Windows, activate the environment with:

```bash
.venv\Scripts\activate
```

Install the required Python packages:

```bash
pip install -r requirements.txt
```

Check that IsoImpact can be started:

```bash
python IsoImpact.py --help
```

## 2. Usage

IsoImpact requires three main inputs:

```text
1. Ensembl transcript IDs for the isoforms to be compared.
2. A matching GTF annotation file.
3. A matching protein FASTA file.
```

The transcript IDs must match the supplied GTF file. The protein FASTA file must contain the protein sequences corresponding to the selected transcripts. For standard human and mouse Ensembl release 110 isoforms, IsoImpact automatically uses the built-in domain-coordinate files in the `data` directory, so users do not need to provide `-d/--domain`.

The GTF file is required because IsoImpact uses it to identify transcript, exon, CDS, UTR, gene, and protein ID annotations for the selected isoforms. The protein FASTA file is required because IsoImpact extracts the corresponding amino acid sequences and calculates protein sequence features from those sequences. For standard human and mouse workflows, the GTF and protein FASTA files should come from Ensembl release 110, because the built-in domain-coordinate files were generated from the same Ensembl release. Using matching files helps ensure that transcript IDs, protein IDs, genomic coordinates, and domain annotations can be linked correctly. These reference GTF and protein FASTA files are large, so they are not included directly in this GitHub repository. Users should download them from Ensembl before running IsoImpact.

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
mkdir -p results/MROH2B

python IsoImpact.py \
  -i ENST00000421030 ENST00000440047 ENST00000413188 ENST00000409996 \
  -g references/Homo_sapiens.GRCh38.110.gtf \
  -f references/Homo_sapiens.GRCh38.pep.all.fa \
  -o results/MROH2B \
  --prefix MROH2B
```

For human transcript IDs beginning with `ENST`, IsoImpact automatically selects:

```text
data/human_domain.csv
```

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
mkdir -p results/mouse_isoforms

python IsoImpact.py \
  -i ENSMUST00000000000 ENSMUST00000000001 \
  -g references/Mus_musculus.GRCm39.110.gtf \
  -f references/Mus_musculus.GRCm39.pep.all.fa \
  -o results/mouse_isoforms \
  --prefix mouse_isoforms
```

Replace the example transcript IDs with the mouse isoforms to be compared. For mouse transcript IDs beginning with `ENSMUST`, IsoImpact automatically selects:

```text
data/mouse_domain.csv
```

## 3. Custom Domain Annotation

This section is only needed for novel isoforms, custom isoforms, non-model organisms, or annotations that are not covered by the built-in human and mouse Ensembl release 110 domain-coordinate files.

For a custom analysis, users need:

```text
1. A CDS-aware GTF file.
2. A matching protein FASTA file.
3. Pfam/PfamScan domain-prediction results for the corresponding protein sequences.
```

The protein IDs in the Pfam/PfamScan result file must match the `protein_id` values in the GTF file. The same protein IDs should also be used in the protein FASTA headers.

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
  -o results/novel \
  --prefix novel
```

Use `-d/--domain` only when providing a custom domain-coordinate file. Standard human and mouse Ensembl release 110 workflows do not require this option.

## 4. Output Files

IsoImpact generates two main output files:

```text
<prefix>_features.csv
<prefix>_figure.png
```

`<prefix>_features.csv` contains transcript metadata, domain summaries, canonical feature values, alternative feature values, and feature-difference values.

`<prefix>_figure.png` contains the domain-mapping visualization and ranked feature-difference plots for the input isoforms.

IsoImpact selects the canonical baseline by parsing `Ensembl_canonical` or `MANE_Select` tags when available. If neither tag is available, it uses the input transcript with the largest genomic span.

## 5. Contact

Please open a GitHub issue for questions, bug reports, or feature requests.
