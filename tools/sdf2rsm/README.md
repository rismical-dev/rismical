# makeparam_rismical.py

A Python tool designed to generate molecular parameter files specifically for **RISMiCal**. It assigns Lennard-Jones parameters (OpenFF Sage) and AM1-BCC partial charges (AmberTools) to a molecule from an SDF file.

## Features
* **RISMiCal Format**: Outputs parameters in a plain-text format compatible with RISMiCal (using `d0` notation for doubles).
* **OpenFF Sage**: Uses the modern Sage force field for high-quality LJ parameters.
* **AmberTools Integration**: Automatically calculates AM1-BCC charges using `sqm`.
* **United Atom (UA) Support**: Optional mode to merge aliphatic hydrogens into their parent carbons.

## Installation
It is recommended to use a Conda-based environment (Miniconda, Anaconda, or Miniforge).

```bash
# 1. Create environment and install dependencies
conda create -n md_param_env -c conda-forge python=3.10 ambertools openff-toolkit openmm rdkit

# 2. Activate the environment
conda activate md_param_env

# 3. Fix compatibility issue (Downgrade setuptools)
conda install -c conda-forge "setuptools<70.0.0"
```

## Usage

### 1. All Atom (AA) Mode
```bash
python makeparam_rismical.py molecule.sdf
```
*Output: `molecule.txt`*

### 2. United Atom (UA) Mode
This mode merges hydrogens bonded to carbons into the parent carbon's charge and removes those hydrogens from the list.
```bash
python makeparam_rismical.py molecule.sdf ua
```
*Output: `molecule.txt`*

## Output Format
The generated text file follows the RISMiCal input style:

**Line 1:** `[Number of Atoms]    [Molecule Name]`  
**Line 2+:** `[Label]  [Sigma]  [Epsilon]  [Charge]  [X]  [Y]  [Z]`

### Units
* **Sigma (σ):** Å (written as `d0`)
* **Epsilon (ε):** J/mol (written as `d0`)
* **Charge (q):** e (written as `d0`)
* **Coordinates (X, Y, Z):** Å

### Example Output
```text
3    mecn
N      3.2500d0   711.280d0  -0.5600d0      1.2608000      0.0000000      0.0000000
C      3.3997d0   457.730d0   0.2000d0     -1.3613000      0.0000000      0.0000000
C      3.3997d0   457.730d0   0.3600d0      0.1006000      0.0000000      0.0000000
```