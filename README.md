# MolecularDynamics
Repo for playing with MD simulations


## Table of Contents

- [MolecularDynamics](#moleculardynamics)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
    - [Prerequisites](#prerequisites)
    - [Gromacs](#gromacs)
  - [Data preparation](#data-preparation)

## Installation

### Prerequisites

### Gromacs
Install GROMACS
```bash
conda install conda-forge::gromacs
```
Confirm installation
```bash
gmx --version
```

install other stuff:
```bash
poetry add biopython
```
## Data preparation
Download a structure from the PDB (e.g. azoreductase from B. subtilis 1NNI), and extract chain A
```bash
# download the structure + cofactors
wget https://files.rcsb.org/download/1NNI.pdb -O data/1NNI.pdb
grep "^ATOM.* A " data/1NNI.pdb > data/1NNI:A.pdb

# download the AlpahFold Structure (no FMN)
wget -O data/AF-O07529-F1-model_v4.pdb https://alphafold.ebi.ac.uk/files/AF-O07529-F1-model_v4.pdb
```
