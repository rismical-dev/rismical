# GMX2RISMiCal

A Python utility to convert GROMACS topology (`.top`) and structure (`.gro`) files into parameter files (`.rsm`) compatible with **RISMiCal**.

## Features

- **Recursive Topology Parsing**: Automatically traces `#include` directives to collect atom types and parameters distributed across multiple files.
- **Automatic Unit Conversion**: 
  - Length: nm to Å (multiplied by 10.0).
  - Energy: kJ/mol to J/mol (multiplied by 1000.0).
- **Flexible LJ Parameter Extraction**: Supports various force field formats (e.g., OPLS-AA, AMBER). It extracts the Lennard-Jones parameters ($\sigma$ and $\epsilon$) from the last columns of the `[ atomtypes ]` section.
- **Robust Processing**: Gracefully skips unneeded include files (such as `posre.itp`) that are not required for parameter extraction, issuing a warning without halting execution.
- **GROMACS Path Integration**: Searches for force field files in the current directory, the path specified by the `GMXDATA` environment variable, and standard Homebrew installation paths.
- **Molecule Centering**: An optional feature to translate the molecule's geometric center to the origin (0, 0, 0). When this feature is used, a `.pdb` file reflecting the translated coordinates is automatically generated.

## Requirements

- Python 3.x
- GROMACS force field files (setup may vary depending on your environment)

## Setup

1. Save the script as `gmx2rismical.py`.
2. (Optional) If you are using standard GROMACS force fields, setting the `GMXDATA` environment variable is recommended so the script can locate them.
   ```bash
   export GMXDATA=/opt/homebrew/share/gromacs
   ```

## Usage

Run the script by providing the GROMACS topology and structure files as arguments. Since the output is printed to the standard output, use redirection (`>`) to create the `.rsm` file.

```bash
python3 gmx2rismical.py topol.top input.gro > molecule.rsm
```

### Using the Centering Option

If you want to align the molecule's geometric center to the origin, use the `-c` or `--centering` flag.

```bash
python3 gmx2rismical.py -c topol.top input.gro > molecule.rsm
```

When this flag is specified:
1. The coordinates inside the generated `.rsm` file will be centered.
2. A PDB file containing the centered coordinates will be automatically generated alongside it (e.g., `input_centering.pdb`).

### Command Line Arguments
- `top_file`: The GROMACS topology file (must contain the `[ atoms ]` section and required `[ atomtypes ]` via include files).
- `gro_file`: The GROMACS structure file corresponding to the topology.
- `-c, --centering` (Optional): Translates the molecule's geometric center to the origin and outputs a `[gro_filename]_centering.pdb` file.

## Output Format

The generated `.rsm` file adheres to the following RISMiCal specifications:

1. **First Line**: `number_of_atoms` and `molecule_name` (derived from the topology file's name).
2. **Subsequent Lines**: One line per atom, outputting the following columns:
   - `atom_label`: Atom name defined in the `[ atoms ]` section.
   - `LJ_sigma`: Lennard-Jones sigma [Å].
   - `LJ_epsilon`: Lennard-Jones epsilon [J/mol].
   - `point_charge`: Partial charge [e].
   - `x`, `y`, `z`: Atomic coordinates [Å].

## Technical Notes

- The script maps atoms under the assumption that the atom order in the `.top` and `.gro` files perfectly matches.
- A `KeyError` will be raised if an atom type defined in `[ atoms ]` cannot be found in any of the included `[ atomtypes ]` sections.
- Warnings regarding skipped include files are printed to standard error (`stderr`), ensuring they do not corrupt or interfere with the redirected `.rsm` file.