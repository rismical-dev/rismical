# GMX2RISMiCal

A Python utility to convert GROMACS topology (`.top`) and structure (`.gro`) files into parameter files compatible with **RISMiCal**.

## Features

- **Recursive Topology Parsing**: Automatically resolves `#include` directives to collect atom types and parameters across multiple files.
- **Unit Conversion**: 
  - Length: nm to Å (multiplied by 10.0).
  - Energy: kJ/mol to J/mol (multiplied by 1000.0).
- **Flexible LJ Parsing**: Handles various force field formats (e.g., OPLS-AA, AMBER) by identifying Lennard-Jones parameters ($\sigma$ and $\epsilon$) from the last columns of the `[ atomtypes ]` section.
- **Robustness**: Gracefully skips missing include files (such as `posre.itp`) that are not required for parameter extraction, preventing unnecessary crashes.
- **GROMACS Integration**: Searches for force field files in the current directory, the path specified by the `GMXDATA` environment variable, and standard Homebrew installation paths.

## Requirements

- Python 3.x
- GROMACS force field files (optional, depending on your system configuration)

## Setup

1. Save the script as `gmx2rismical.py`.
2. (Optional) Set the `GMXDATA` environment variable to point to your GROMACS data directory to ensure the script can find standard force fields.
   ```bash
   export GMXDATA=/opt/homebrew/share/gromacs
   ```

## Usage

Run the script by providing the GROMACS topology file and the coordinate file as arguments. You can redirect the output to a `.prm` file.

```bash
python3 gmx2rismical.py topol.top protein.gro > molecule.prm
```

### Command Line Arguments
- `topol.top`: The GROMACS topology file (should contain `[ atoms ]` and include files for `[ atomtypes ]`).
- `protein.gro`: The GROMACS coordinate file matching the topology.

## Output Format

The output follows the RISMiCal parameter file specification:

1. **First Line**: `number_of_atoms` and `molecule_name` (extracted from the `.top` filename).
2. **Subsequent Lines**: One line per atom with the following columns:
   - `atom_label`: Atom name from the `[ atoms ]` section.
   - `LJ_sigma`: Lennard-Jones sigma in Angstroms [Å].
   - `LJ_epsilon`: Lennard-Jones epsilon in Joules per mole [J/mol].
   - `point_charge`: Atomic partial charge in elementary charge units [e].
   - `x`, `y`, `z`: Atomic coordinates in Angstroms [Å].

## Technical Notes

- The script maps atoms between the `.top` and `.gro` files based on their order. Ensure the files correspond to the same system.
- If an atom type defined in `[ atoms ]` is not found in any of the included `[ atomtypes ]` sections, the script will raise a `KeyError`.
- Warnings about skipped include files are printed to `stderr` and will not affect the redirected output file.

