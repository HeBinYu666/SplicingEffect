# IsoImpact

IsoImpact is a command-line pipeline that automates the functional evaluation of protein-coding splice isoforms. Starting from isoform IDs, a GTF annotation file, a proteome FASTA file, and a matching domain-coordinate CSV file, IsoImpact integrates domain mapping, multi-isoform domain visualization and comparison, protein feature extraction, and feature-difference analysis into a unified automated workflow.

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

The transcript IDs (`-i`) must match the supplied GTF file. The GTF file (`-g`) provides transcript, exon, CDS, UTR, gene, and protein ID annotations for the selected isoforms. The protein FASTA file (`-f`) provides the corresponding amino acid sequences. The domain-coordinate CSV file (`-d`) provides the protein-domain annotations used for domain comparison and visualization.

IsoImpact includes pre-built human and mouse domain-coordinate CSV files under the `data/` directory. These files were generated from Ensembl release 110, so the human and mouse examples below use Ensembl release 110 GTF and protein FASTA files to keep transcript IDs, protein IDs, genomic coordinates, and domain annotations consistent.

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

Because this example compares human isoforms, it uses the human domain-coordinate file `data/human_domain.csv`.

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

`IsoImpact_features.csv` contains gene, transcript, protein, coding biotype, protein feature, genomic span, and domain comparison results. An example feature-matrix preview from four human MROH7 isoforms is shown below. The complete example table is available as [`docs/example_output/IsoImpact_features.csv`](docs/example_output/IsoImpact_features.csv).

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

Replace the example transcript IDs with the mouse isoforms to be compared. Because this example compares mouse isoforms, it uses the mouse domain-coordinate file `data/mouse_domain.csv`.

## 3. Building Domain Annotation Files for Novel Isoforms

For novel isoforms or other isoforms that are not included in the provided domain-coordinate CSV files, users can use `scripts/build_custom_domain.R` to build a matching domain-coordinate CSV file for IsoImpact. Here, novel isoforms include newly identified human or mouse isoforms and isoforms from other species.

For this process, users need:

```text
1. A GTF file that contains the target isoforms and their CDS records, including transcript_id and protein_id annotations.
2. Pfam/PfamScan domain-prediction results for the protein sequences of the same target isoforms.
```

The GTF file must contain CDS records because the helper script maps protein-domain intervals back to genomic CDS coordinates. The protein IDs in the Pfam/PfamScan result file must match the `protein_id` values in the GTF file; otherwise, the script cannot link the predicted domains to the corresponding isoforms.

Install the R packages required by the helper script:

```r
install.packages("BiocManager")
BiocManager::install(c("ensembldb", "GenomicRanges"))
install.packages("dplyr")
```

Build an IsoImpact-compatible domain-coordinate file:

```bash
Rscript scripts/build_custom_domain.R \
  --gtf your_novel_isoforms.gtf \
  --pfam your_novel_pfam_results.txt \
  --out novel_custom_domain.csv
```

Then run IsoImpact with the generated domain-coordinate file:

```bash
python IsoImpact.py \
  -i novel_tx_1 novel_tx_2 \
  -g your_novel_isoforms.gtf \
  -f your_novel_proteins.fa \
  -d novel_custom_domain.csv \
  -o results/novel
```

## 4. Contact

If you have questions about installing or running IsoImpact, or if you find an error in the code or documentation, please contact:

- Bin-Yu He: hzb022119@163.com
- Hong-Dong Li: hongdong@csu.edu.cn

You can also click the **Issues** tab at the top of this GitHub page to report problems.
