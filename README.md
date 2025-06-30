# Ursina Proteins

![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2FHarrisonTCodes%2Fursina-proteins%2Frefs%2Fheads%2Fmain%2Fpyproject.toml)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
![Pre-Commit Check](https://github.com/HarrisonTCodes/ursina-proteins/actions/workflows/pre-commit.yaml/badge.svg)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://pre-commit.com/)
[![Conventional Commits](https://img.shields.io/badge/Conventional%20Commits-1.0.0-%23FE5196?logo=conventionalcommits&logoColor=white)](https://conventionalcommits.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A Python package for rendering protein structures of PDB format in 3D using [Ursina](https://www.ursinaengine.org/), as featured in [Ursina's sample projects](https://www.ursinaengine.org/samples.html).

![Example proteins](./assets/example.png)

## Installation
The package is published on [PyPI](https://pypi.org/project/ursina-proteins/), and can be installed with `pip install ursina-proteins`.
You can also clone the repo down and install dependencies with [uv](https://docs.astral.sh/uv/).
```bash
# Clone the repo
git clone https://github.com/HarrisonTCodes/ursina-proteins.git
cd ursina_proteins

# Install dependencies with uv
uv sync
```

## Usage
You can use the library in an existing Ursina project by importing the `Protein` class and creating an instance from a PDB file. You can render any proteins in the PDB file format, and can find many available for download at the [RCSB Protein Data Bank](https://www.rcsb.org/).
```python
from ursina_proteins.protein import Protein

Protein("/path/to/file.pdb")
```
You can also test the library out by running [demo.py](https://github.com/HarrisonTCodes/ursina-proteins/blob/main/src/demo.py). This script renders a simple scene with an example protein ([insulin](https://www.rcsb.org/structure/3I40)).
```bash
uv run src/demo.py
```

## Contributions
Contributions are welcome. Please enable pre-commit hooks to catch and fix formatting/linting issues locally before raising a PR.
```bash
# Enable pre-commit hooks
uv run pre-commit install
```
