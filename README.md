# Ursina Proteins

![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2FHarrisonTCodes%2Fursina-proteins%2Frefs%2Fheads%2Fmain%2Fpyproject.toml)
[![Poetry](https://img.shields.io/endpoint?url=https://python-poetry.org/badge/v0.json)](https://python-poetry.org/)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://pre-commit.com/)
[![black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Conventional Commits](https://img.shields.io/badge/Conventional%20Commits-1.0.0-%23FE5196?logo=conventionalcommits&logoColor=white)](https://conventionalcommits.org)

A Python package for visualizing protein structures of PDB format in 3D using [Ursina](https://www.ursinaengine.org/).

## Installation
To install the library locally, clone the repo down and install dependencies with Poetry.
```bash
# Clone the repo
git clone https://github.com/HarrisonTCodes/ursina-proteins.git
cd ursina_proteins

# Install with poetry
poetry install
```

## Usage
You can use the library in an existing Ursina project by importing the Protein class and creating an instance from a PDB file.
```python
from ursina_proteins.protein import Protein

Protein("/path/to/file.pdb")
```
You can also test the library out and render a protein by running the `demo.py`. You can render any proteins in the PDB file format, such as [hemoglobin](https://doi.org/10.2210/pdb1a3n/pdb) or this [small example file](https://gist.github.com/cstein/6699200).
```bash
poetry run python src/demo.py
# You will then be prompted to enter the path to a PDB file
```

## Contributions
Contributions are welcome. Please use conventional commits and enable pre-commit hooks.
```bash
# Enabled pre-commit hooks
poetry run pre-commit install
```
