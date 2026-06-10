# Generate_LoQI

This folder contains the local wrapper used to generate ligand conformers with [LoQI](https://github.com/isayevlab/LoQI), an external package for low-energy QM-informed conformer generation.

## Files

### `generate_loqi_generic.py`

General batch-processing script for running LoQI on ligand SMILES CSV files.

Functionality:

* Reads one or more CSV files containing ligand IDs and SMILES strings.
* Runs LoQI conformer generation for each molecule.
* Saves generated conformers as `.sdf` files.
* Saves per-molecule metadata as `.json` files.
* Writes batch logs, result summaries, failed-molecule records, and checkpoints.
* Supports resuming interrupted runs from checkpoint files.
* Can optionally convert generated `.sdf` conformers to `.xyz`.

Expected CSV columns:

```text
id, smiles
```

Example usage:

```bash
python Generate_LoQI/generate_loqi_generic.py \
  --csvs input/ML-CN-Smiles.csv input/ML-NN-Smiles.csv \
  --output-base Generate_LoQI/output \
  --n-conformers 12 \
  --batch-size 50
```

### `generate_loqi.ipynb`

Notebook interface for running `generate_loqi_generic.py` with project-specific paths. This notebook is mainly used to configure local paths, launch LoQI generation, and inspect outputs during development.

## External dependency: LoQI

This code relies on the external LoQI package. The LoQI source code, checkpoint files, and required data files are not authored in this project.

Recommended local structure:

```text
Generate_LoQI/
├── generate_loqi_generic.py
├── generate_loqi.ipynb
├── requirements-loqi-mac.txt
├── LoQI/
└── output/
```

The `LoQI/` folder should contain the cloned LoQI repository.

Required data/checkpoint files must be downloaded from https://doi.org/10.1184/R1/31441570, and placed with the following layout:

```text
LoQI/
└── data/
    ├── loqi.ckpt
    ├── loqi_flow.ckpt
    └── chembl3d_stereo/
        ├── processed/
        ...
```

## Outputs

Each successful molecule produces:

```text
.sdf   generated 3D conformer structures
.json  metadata for the LoQI run
```

Batch-level outputs include logs, checkpoint files, summary CSVs, and failed-molecule CSVs.

## Citations

This project uses LoQI and its underlying Megalodon architecture to generate molecular conformers for predictive modeling.

```bibtex
@article{nikitin2025scalable,
  title={Scalable Low-Energy Molecular Conformer Generation with Quantum Mechanical Accuracy},
  author={Nikitin, Filipp and Anstine, Dylan M and Zubatyuk, Roman and Paliwal, Saee Gopal and Isayev, Olexandr},
  year={2025}
}

@article{reidenbach2025applications,
  title={Applications of Modular Co-Design for De Novo 3D Molecule Generation},
  author={Reidenbach, Danny and Nikitin, Filipp and Isayev, Olexandr and Paliwal, Saee},
  journal={arXiv preprint arXiv:2505.18392},
  year={2025}
}
```
