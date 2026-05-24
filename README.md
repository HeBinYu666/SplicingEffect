# IsoImpact

IsoImpact is an automated pipeline to evaluate the functional impact of alternative splicing based on isoform sequences. It integrates multi-isoform domain visualization and comparison, protein feature extraction, and feature-difference analysis into a unified automated workflow.

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

## 2. Command-line options

IsoImpact uses the following command-line inputs:

```text
-i / --isoform   Ensembl transcript IDs for the isoforms to be compared. At least two transcript IDs are required.
-g / --gtf       Path to the GTF annotation file.
-f / --fasta     Path to the protein FASTA file.
-d / --domain    Path to the domain-coordinate CSV file.
-o / --outdir    Directory where output files will be saved. The default is the current directory.
```

The transcript IDs (`-i`) must match the supplied GTF file. The GTF file (`-g`) provides transcript, exon, CDS, UTR, gene, and protein ID annotations for the selected isoforms. The protein FASTA file (`-f`) provides the corresponding amino acid sequences. The domain-coordinate CSV file (`-d`) provides the protein-domain annotations used for domain comparison and visualization. For standard human and mouse analyses, IsoImpact provides pre-built domain-coordinate CSV files under the `data/` directory.

These pre-built human and mouse domain-coordinate CSV files were generated from Ensembl release 110, so the human and mouse examples below use Ensembl release 110 GTF and protein FASTA files to keep transcript IDs, protein IDs, genomic coordinates, and domain annotations consistent. Other Ensembl releases can also be used after generating a matching domain-coordinate CSV file from the same annotation version.

### 2.1 Example usage using human isoforms

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

Step 2. Run IsoImpact with four human Ensembl transcript IDs from **MROH7** (**ENSG00000184313**):

```bash
python IsoImpact.py \
  -i ENST00000421030 ENST00000440047 ENST00000413188 ENST00000409996 \
  -g references/Homo_sapiens.GRCh38.110.gtf \
  -f references/Homo_sapiens.GRCh38.pep.all.fa \
  -d data/human_domain.csv \
  -o results/human_isoforms
```

Because this example compares human isoforms, it uses the human domain-coordinate file `data/human_domain.csv` by default.

After the command finishes, IsoImpact writes the output figure to:

```text
results/human_isoforms/IsoImpact_figure.png
```

`IsoImpact_figure.png` contains the domain-mapping visualization and ranked feature-difference plot for the input isoforms. An example output figure from the four human MROH7 isoforms is shown below. The example figure is available as [`docs/example_output/IsoImpact_figure.png`](docs/example_output/IsoImpact_figure.png).

![IsoImpact example output figure](docs/example_output/IsoImpact_figure.png)

IsoImpact also writes the feature matrix to:

```text
results/human_isoforms/IsoImpact_features.csv
```

`IsoImpact_features.csv` contains gene and transcript information, protein features, genomic span, and domain comparison results. An example feature-matrix preview from four human MROH7 isoforms is shown below. The complete example table is available as [`docs/example_output/IsoImpact_features.csv`](docs/example_output/IsoImpact_features.csv).

| Gene ID | Canonical transcript | Alternative transcript | Alternative domains | Lost domains | Gained domains | Total domain change |
| --- | --- | --- | --- | --- | --- | ---: |
| ENSG00000184313 | ENST00000421030 | ENST00000440047 | HEAT_Maestro_2; HEAT_Maestro | HEAT_MROH2B_C | None | 1 |
| ENSG00000184313 | ENST00000421030 | ENST00000413188 | HEAT_Maestro_2 | HEAT_Maestro; HEAT_MROH2B_C | None | 2 |
| ENSG00000184313 | ENST00000421030 | ENST00000409996 | HEAT_Maestro_2; HEAT_Maestro; HEAT_MROH2B_C | None | None | 0 |

### 2.2 Example usage using mouse isoforms

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

Replace the example transcript IDs with the mouse isoforms to be compared. Because this example compares mouse isoforms, it uses the mouse domain-coordinate file `data/mouse_domain.csv` by default.

## 3. Building Domain-Coordinate CSV Files for Other Ensembl Releases

For isoforms annotated in Ensembl releases other than release 110, users can use the `build_custom_domain.R` script in the `scripts/` directory to build a matching domain-coordinate CSV file. The generated CSV file should be used together with the GTF and protein FASTA files from the same Ensembl release.

For this process, users need:

```text
1. A GTF file from the target Ensembl release.
2. A protein FASTA file from the same Ensembl release.
3. PfamScan domain-prediction results for the corresponding protein sequences.
```

The GTF file is used by IsoImpact for transcript and CDS coordinates. The PfamScan result file should retain the original output structure generated by PfamScan, without manual modification.

Install the R packages required by the helper script:

```r
install.packages("BiocManager")
BiocManager::install(c("AnnotationHub", "ensembldb", "GenomicRanges"))
install.packages("dplyr")
```

The example below uses two NOTCH2 isoforms from Ensembl release 109. To avoid running PfamScan on the full human proteome, this repository provides a small protein FASTA subset from Ensembl release 109:

```text
docs/example_release109/Homo_sapiens.GRCh38.pep.release109.example.fa
```

Step 1. Download the Ensembl release 109 human GTF file. To save time, the example workflow below uses the small protein FASTA subset already provided in this repository. Users only need to download the complete Ensembl release 109 protein FASTA file if they want to run PfamScan on additional protein sequences from the same release.

```bash
mkdir -p references

curl -L -o references/Homo_sapiens.GRCh38.109.gtf.gz \
  https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

gunzip references/Homo_sapiens.GRCh38.109.gtf.gz
```

Optional command for downloading the complete Ensembl release 109 protein FASTA file:

```bash
curl -L -o references/Homo_sapiens.GRCh38.pep.all.release109.fa.gz \
  https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz

gunzip references/Homo_sapiens.GRCh38.pep.all.release109.fa.gz
```

Step 2. Install PfamScan, prepare the Pfam database, and run PfamScan on the provided small protein FASTA subset. The Pfam database is large, so it is not included in this repository. If PfamScan and the Pfam database are already available on your computer or server, skip the installation commands and use the existing paths when running PfamScan.

Install PfamScan and prepare the Pfam database:

```bash
conda create -n pfamscan -c conda-forge -c bioconda pfam_scan hmmer -y
conda activate pfamscan

mkdir -p references/pfam

curl -L -o references/pfam/Pfam-A.hmm.gz \
  https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

curl -L -o references/pfam/Pfam-A.hmm.dat.gz \
  https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz

curl -L -o references/pfam/active_site.dat.gz \
  https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz

gunzip -f references/pfam/Pfam-A.hmm.gz
gunzip -f references/pfam/Pfam-A.hmm.dat.gz
gunzip -f references/pfam/active_site.dat.gz

hmmpress references/pfam/Pfam-A.hmm
```

Run PfamScan:

```bash
mkdir -p results/release109_notch2

pfam_scan.pl \
  -fasta docs/example_release109/Homo_sapiens.GRCh38.pep.release109.example.fa \
  -dir references/pfam \
  -outfile results/release109_notch2/NOTCH2_release109_pfam_results.txt
```

If PfamScan and the Pfam database are already installed on your computer, replace the program and database paths:

```bash
perl /path/to/pfam_scan.pl \
  -fasta docs/example_release109/Homo_sapiens.GRCh38.pep.release109.example.fa \
  -dir /path/to/pfam_database \
  -outfile results/release109_notch2/NOTCH2_release109_pfam_results.txt
```

Step 3. Build an IsoImpact-compatible domain-coordinate file:

```bash
Rscript scripts/build_custom_domain.R \
  --gtf references/Homo_sapiens.GRCh38.109.gtf \
  --pfam results/release109_notch2/NOTCH2_release109_pfam_results.txt \
  --out results/release109_notch2/NOTCH2_release109_domain.csv \
  --ensembl 109
```

Step 4. Run IsoImpact with the generated domain-coordinate file:

```bash
python IsoImpact.py \
  -i ENST00000256646 ENST00000652302 \
  -g references/Homo_sapiens.GRCh38.109.gtf \
  -f docs/example_release109/Homo_sapiens.GRCh38.pep.release109.example.fa \
  -d results/release109_notch2/NOTCH2_release109_domain.csv \
  -o results/release109_notch2
```

After the command finishes, IsoImpact writes the output figure to:

```text
results/release109_notch2/IsoImpact_figure.png
```

`IsoImpact_figure.png` contains the domain-mapping visualization and ranked feature-difference plot for the input isoforms.

[Add the generated `IsoImpact_figure.png` here.]

IsoImpact also writes the feature matrix to:

```text
results/release109_notch2/IsoImpact_features.csv
```

`IsoImpact_features.csv` contains gene and transcript information, protein features, genomic span, and domain comparison results.

[Add a preview of the generated `IsoImpact_features.csv` here.]

## 4. Contact

If you have questions about installing or running IsoImpact, or if you find an error in the code or documentation, please contact:

- Bin-Yu He: hzb022119@163.com
- Hong-Dong Li: hongdong@csu.edu.cn

You can also click the **Issues** tab at the top of this GitHub page to report problems.
