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

## 2. Usage

IsoImpact uses the following command-line inputs:

```text
-i / --isoform   Ensembl transcript IDs for the isoforms to be compared. At least two transcript IDs are required.
-g / --gtf       Path to the matching GTF annotation file.
-f / --fasta     Path to the matching protein FASTA file.
-d / --domain    Path to the matching domain-coordinate CSV file.
-o / --outdir    Directory where output files will be saved. The default is the current directory.
```

The transcript IDs must match the supplied GTF file. The GTF file provides transcript, exon, CDS, UTR, gene, and protein ID annotations for the selected isoforms. The protein FASTA file provides the corresponding amino acid sequences. The domain-coordinate CSV file provides the protein-domain annotations used for domain comparison and visualization.

The provided human and mouse domain-coordinate CSV files were built from Ensembl release 110. Therefore, the human and mouse examples below use Ensembl release 110 GTF and protein FASTA files so that transcript IDs, protein IDs, genomic coordinates, and domain annotations match each other. Users can use another annotation version by building a matching domain-coordinate CSV file with the helper script described in Section 3.

### 2.1 Human

Step 1. Download the Ensembl release 110 human GTF and protein FASTA files:

```bash
mkdir -p references

curl -L -o references/Homo_sapiens.GRCh38.110.gtf.gz \
  https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

curl -L -o references/Homo_sapiens.GRCh38.pep.all.fa.gz \
  https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz

gunzip references/Homo_sapiens.GRCh38.110.gtf.gz
gunzip references/Homo_sapiens.GRCh38.pep.all.fa.gz
```

Step 2. Run IsoImpact with human Ensembl transcript IDs:

```bash
python IsoImpact.py \
  -i ENST00000421030 ENST00000440047 ENST00000413188 ENST00000409996 \
  -g references/Homo_sapiens.GRCh38.110.gtf \
  -f references/Homo_sapiens.GRCh38.pep.all.fa \
  -d data/human_domain.csv \
  -o results/human_isoforms
```

This example uses `data/human_domain.csv`, which was built from Ensembl release 110.

After the command finishes, IsoImpact writes two main output files to the output directory:

```text
results/human_isoforms/IsoImpact_features.csv
results/human_isoforms/IsoImpact_figure.png
```

`IsoImpact_features.csv` contains gene, transcript, protein, coding biotype, protein feature, genomic span, and domain comparison results. An example feature-table preview from four human MROH2B isoforms is shown below. The complete example table is available as [`docs/example_output/IsoImpact_features.csv`](docs/example_output/IsoImpact_features.csv).

| Alternative transcript | Alternative domains | Lost domains | Gained domains | Total domain change |
| --- | --- | --- | --- | ---: |
| ENST00000440047 | HEAT_Maestro_2; HEAT_Maestro | HEAT_MROH2B_C | None | 1 |
| ENST00000413188 | HEAT_Maestro_2 | HEAT_Maestro; HEAT_MROH2B_C | None | 2 |
| ENST00000409996 | HEAT_Maestro_2; HEAT_Maestro; HEAT_MROH2B_C | None | None | 0 |

`IsoImpact_figure.png` contains the domain-mapping visualization and ranked feature-difference plot for the input isoforms. An example output figure from the same four human MROH2B isoforms is shown below. The example figure is available as [`docs/example_output/IsoImpact_figure.png`](docs/example_output/IsoImpact_figure.png).

![IsoImpact example output figure](docs/example_output/IsoImpact_figure.png)

### 2.2 Mouse

Step 1. Download the Ensembl release 110 mouse GTF and protein FASTA files:

```bash
mkdir -p references

curl -L -o references/Mus_musculus.GRCm39.110.gtf.gz \
  https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz

curl -L -o references/Mus_musculus.GRCm39.pep.all.fa.gz \
  https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz

gunzip references/Mus_musculus.GRCm39.110.gtf.gz
gunzip references/Mus_musculus.GRCm39.pep.all.fa.gz
```

Step 2. Run IsoImpact with mouse Ensembl transcript IDs:

```bash
python IsoImpact.py \
  -i ENSMUST00000112701 ENSMUST00000134301 \
  -g references/Mus_musculus.GRCm39.110.gtf \
  -f references/Mus_musculus.GRCm39.pep.all.fa \
  -d data/mouse_domain.csv \
  -o results/mouse_isoforms
```

Replace the example transcript IDs with the mouse isoforms to be compared. This example uses `data/mouse_domain.csv`, which was built from Ensembl release 110.

## 3. Building Domain Annotation Files for New Isoforms

The human and mouse examples above use the provided domain-coordinate CSV files, so no R package is required for those analyses. The R script in this repository is provided for isoforms that are not recorded in the provided domain-coordinate files. In this situation, users need to prepare their own annotation files and domain-prediction results, then use the R script to build a matching domain-coordinate CSV file before running IsoImpact.

Here, new isoforms means isoforms that are not covered by the domain-coordinate CSV file used for the analysis. If the isoforms are not recorded in the GTF and protein FASTA files, users also need to provide matching GTF and protein FASTA files for those isoforms. To build the domain-coordinate CSV file, the helper script uses the GTF file and Pfam/PfamScan results.

For this process, users need:

```text
1. A CDS-aware GTF file.
2. Pfam/PfamScan domain-prediction results for the corresponding protein sequences.
```

The GTF file must contain CDS records, because the helper script maps protein-domain intervals back to genomic CDS coordinates. The protein IDs in the Pfam/PfamScan result file must match the `protein_id` values in the GTF file.

Install the R packages required by the helper script:

```r
install.packages("BiocManager")
BiocManager::install(c("ensembldb", "GenomicRanges"))
install.packages("dplyr")
```

Build an IsoImpact-compatible domain-coordinate file:

```bash
Rscript scripts/build_custom_domain.R \
  --gtf your_isoforms.gtf \
  --pfam your_pfam_results.txt \
  --out custom_domain.csv
```

Then run IsoImpact with the generated domain-coordinate file:

```bash
python IsoImpact.py \
  -i novel_tx_1 novel_tx_2 \
  -g your_isoforms.gtf \
  -f your_proteins.fa \
  -d custom_domain.csv \
  -o results/novel
```

## 4. Contact

Name:

Email:

Name:

Email:
