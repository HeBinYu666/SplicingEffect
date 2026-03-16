# SplicingEffect 🧬

**SplicingEffect** is a fully automated computational pipeline for evaluating and visualizing the functional consequences of alternative splicing events at the isoform level. 

## Features
- **Single Command Execution**: Automates canonical identification, domain mapping, and feature extraction without tedious manual operations.
- **Dynamic Visualization**: Generates a high-quality, multi-track mapping of protein domains onto physical genomic coordinates.
- **Ultra-high Dimensional Features**: Extracts nearly 3,000 features per transcript pair, including physicochemical properties, 147-dim CTD patterns, and spatial PseAAC/Autocorrelation features for downstream machine learning.

## Quick Start (with Example Data)

We provide a mini-dataset in the `example_data/` directory for quick testing. You can run the entire pipeline with a single command:

```bash
python SplicingEffect.py -i ENST00000375213 ENST00000338983 ENST00000409176
