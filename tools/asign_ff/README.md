# assign_ff

**Interactive force-field parameter generator for small molecules**, using
[AmberTools](https://ambermd.org/AmberTools.php) (`antechamber`, `parmchk2`,
`tleap`) and optionally [Gaussian 16](https://gaussian.com/) as external
programs.

Given an arbitrary small molecule, `assign_ff.py` assigns GAFF or GAFF2 atom
types, computes partial charges (AM1-BCC or RESP), and writes parameter /
topology / coordinate files for **Amber**, **Gromacs**, and the **3D-RISM**
`.rsm` format.

---

## Features

- **Multiple input formats**, auto-detected from extension and content:
  `pdb`, `xyz`, `mol2`, `sdf` (`.mol`), `gjf` (`.com`, `.gau`),
  `gout` (`.log`, `.out`).
- **Two force fields**: GAFF or GAFF2.
- **Two charge methods**:
  - **AM1-BCC** &mdash; fast, no QM required.
  - **RESP** &mdash; HF/6-31G\* ESP via Gaussian 16, followed by RESP fit.
    If the input is a Gaussian output that already contains ESP fit points,
    the QM step is skipped and RESP is read directly.
- **Multiple output formats**, selectable in any combination:
  - **Amber**: `<stem>.prep` + `<stem>.frcmod`
  - **Gromacs**: `<stem>_GMX.top` + `<stem>_GMX.gro` (via `tleap` + ParmEd)
  - **3D-RISM**: `<stem>.rsm` (built directly from `parmchk2 -a Y` + the
    charged mol2; no `tleap` / ParmEd dependency)
- **Geometry control**:
  - Gaussian is invoked with `NoSymm`, so the input frame is preserved.
  - You can optionally translate the center of mass to the origin
    (pure translation, no rotation).
- **Configurable residue name** (up to 3 alphanumeric characters,
  default `MOL`).

---

## Requirements

| Component        | When needed                                | Notes                                                                  |
| ---------------- | ------------------------------------------ | ---------------------------------------------------------------------- |
| Python           | always                                     | &ge; 3.7                                                               |
| `antechamber`    | always                                     | from AmberTools                                                        |
| `parmchk2`       | always                                     | from AmberTools                                                        |
| `tleap`          | only for Gromacs output                    | from AmberTools                                                        |
| ParmEd           | only for Gromacs output                    | bundled with AmberTools; otherwise `pip install parmed`                |
| Gaussian 16 (`g16`) | only for RESP, and only when QM is needed | not needed if the input is a `gout` that already contains ESP fit data |

All AmberTools executables must be on `PATH` (`source ${AMBERHOME}/amber.sh`).
`g16` must be set up per its own documentation (`g16root`, profile script).

The script has **no Python package dependencies** beyond the standard library
unless you request Gromacs output (which needs ParmEd).

---

## Installation

```bash
git clone https://github.com/<you>/assign_ff.git
cd assign_ff
chmod +x assign_ff.py            # optional
```

There is nothing to build &mdash; `assign_ff.py` is a single-file script.

---

## Quick start

```bash
python assign_ff.py molecule.mol2
```

If no path is given on the command line, the script asks for one
interactively. You will then be prompted, in order, for:

1. **Residue name** (default `MOL`)
2. **Total charge** (integer)
3. **Spin multiplicity** (integer, e.g. `1` for singlet)
4. **Force field**: GAFF or GAFF2
5. **Charge method**: AM1-BCC or RESP
   - For RESP, if a QM job is required, you are also asked whether to run
     a geometry optimization (default: **no**).
6. **Output format(s)**: any combination of Amber, Gromacs, 3D-RISM
7. **Geometry frame**: keep original input geometry, or translate the mass
   center to the origin

Outputs are written to `./<stem>_param/`, where `<stem>` is the input
filename without its extension.

---

## Example session

```text
$ python assign_ff.py glycine.gjf

Input file     : /path/to/glycine.gjf
Detected format: gjf

Enter residue name (alphanumeric, up to 3 chars, default MOL): GLY

Enter total molecular charge (integer): 0
Enter spin multiplicity (integer, e.g. 1 = singlet, 2 = doublet): 1

Select the force field to use:
    [1] GAFF  (General Amber Force Field)
    [2] GAFF2 (General Amber Force Field v2)
    Enter a number: 2

Select the charge method:
    [1] AM1-BCC (fast, no QM required)
    [2] RESP    (HF/6-31G* ESP, uses Gaussian 16)
    Enter a number: 2

Select output format(s) (multiple allowed):
    [1] Amber   (prep + frcmod)
    [2] Gromacs (top + gro)
    [3] 3D-RISM (.rsm)
    Enter number(s) separated by commas (e.g. 1,3): 1,3

Select the geometry for the output coordinates:
    [1] Keep the original input geometry
    [2] Translate the mass center to the origin
    Enter a number: 1
...
=== Done ===
Residue name    : GLY
Output directory: /path/to/glycine_param
  [OK ] glycine.prep
  [OK ] glycine.frcmod
  [OK ] glycine.rsm
```

---

## Input formats

The format is detected from the extension and confirmed by inspecting the
file contents.

| Format | Typical extensions      | Notes                                                                                          |
| ------ | ----------------------- | ---------------------------------------------------------------------------------------------- |
| PDB    | `.pdb`, `.ent`          | Bonds may be perceived from distances.                                                         |
| XYZ    | `.xyz`                  | `antechamber` cannot read XYZ directly, so the script converts it to a minimal PDB internally. |
| MOL2   | `.mol2`                 | TRIPOS MOL2 with explicit bonds (recommended for complex molecules).                           |
| SDF    | `.sdf`, `.mol`          | MDL molfile (V2000 / V3000).                                                                   |
| GJF    | `.gjf`, `.com`, `.gau`  | Gaussian input.                                                                                |
| GOUT   | `.gout`, `.log`, `.out` | Gaussian output. If it contains ESP fit data (`iop(6/33=2)`), RESP is read directly.           |

---

## RESP workflow

When RESP is selected:

1. If the input is a `gout` whose log contains `ESP Fit Center` entries
   (i.e. it was produced with `iop(6/33=2)`), the script reads RESP charges
   from it directly. **No QM job is run.**
2. Otherwise the script generates a Gaussian input via `antechamber -fo gcrt`
   with the route line:
   ```
   #P HF/6-31G* [Opt] SCF=Tight NoSymm Pop=MK IOp(6/33=2,6/41=10,6/42=17)
   ```
   `Opt` is included only if you answered "yes" to the optimization prompt
   (default: no). `NoSymm` is **always** included to prevent Gaussian from
   reorienting the molecule to its standard orientation, so that the output
   coordinates match the input frame.
3. `g16` is then invoked, and the resulting `gout` is fed to
   `antechamber -c resp`.

You can change the QM level, memory, and parallelism by editing the
constants at the top of `assign_ff.py`:

```python
QM_METHOD = "HF/6-31G*"
QM_EXTRA  = "SCF=Tight NoSymm Pop=MK IOp(6/33=2,6/41=10,6/42=17)"
G16_NPROC = 4
G16_MEM   = "2GB"
```

---

## Output: 3D-RISM `.rsm` format

```
 $UDATA
       NATOM  LABEL
ATOM1   SIGMA[A]   EPSILON[J/mol]   CHARGE[e]   X[A]   Y[A]   Z[A]
...
 $END
```

- Line 1: `" $UDATA"` (one leading space).
- Line 2: atom count as a 10-digit integer, followed by the label (the
  input filename without its extension).
- Subsequent lines: atom label, Lennard-Jones diameter sigma in angstrom,
  epsilon in J/mol, partial charge in e, then Cartesian coordinates in
  angstrom.
- Last line: `" $END"` (one leading space).

sigma is converted from Amber's `Rmin/2` as `sigma = (Rmin/2) * 2^(5/6)`.
epsilon is converted from kcal/mol to J/mol with the factor 4184.

Lennard-Jones parameters for every atom type that appears in the molecule
are obtained from `parmchk2 -a Y`, which writes the complete LJ table to a
single frcmod (including pre-existing types from `gaff.dat` / `gaff2.dat`).
Atom names, charges, and coordinates are taken from the charged mol2
produced by `antechamber`. **This path does not require `tleap` or ParmEd.**

---

## Output: Amber and Gromacs

- **Amber**: `<stem>.prep` (prepi format, residue name as configured) and
  `<stem>.frcmod` (missing parameters only, from `parmchk2`).
- **Gromacs**: a complete prmtop is first built with `tleap`:
  ```
  source leaprc.gaff[2]
  <RES> = loadmol2 <stem>.mol2
  loadamberparams <stem>.frcmod
  saveamberparm <RES> <stem>.prmtop <stem>.inpcrd
  ```
  then converted with ParmEd to `<stem>_GMX.top` and `<stem>_GMX.gro`.

---

## Geometry frame

Two factors can change the coordinate frame between the input file and the
written outputs:

1. **Gaussian's reorientation**: by default Gaussian transforms the molecule
   into its standard orientation. The script disables this with `NoSymm`,
   so the input frame is preserved.
2. **Optional center-of-mass translation**: if you select
   "translate the mass center to the origin", the charged mol2 is translated
   so that its mass-weighted center sits at `(0, 0, 0)`. This is a pure
   translation; no rotation is applied, and bonds and charges are not
   modified. The translation is applied at the mol2 stage, so all subsequent
   outputs (Amber prep, Gromacs gro, RISM rsm) share the same frame.

---

## Files written

For a run with input `<stem>.<ext>`, the output directory `<stem>_param/`
will contain a subset of the following, depending on the options chosen:

| File                  | Always | Amber | Gromacs | RISM | RESP w/ QM |
| --------------------- | :----: | :---: | :-----: | :--: | :--------: |
| `<stem>.mol2`         |   x    |       |         |      |            |
| `<stem>_resp.gjf`     |        |       |         |      |     x      |
| `<stem>_resp.gout`    |        |       |         |      |     x      |
| `<stem>.prep`         |        |   x   |         |      |            |
| `<stem>.frcmod`       |        |   x   |    x    |      |            |
| `<stem>.prmtop`       |        |       |    x    |      |            |
| `<stem>.inpcrd`       |        |       |    x    |      |            |
| `<stem>_GMX.top`      |        |       |    x    |      |            |
| `<stem>_GMX.gro`      |        |       |    x    |      |            |
| `<stem>_all.frcmod`   |        |       |         |  x   |            |
| `<stem>.rsm`          |        |       |         |  x   |            |

`antechamber` also leaves auxiliary files (`ANTECHAMBER_*`, `ATOMTYPE.INF`,
`punch`, `qout`, `esout`, etc.) in the working directory.

---

## Troubleshooting

**`antechamber` / `parmchk2` / `tleap` not found.**
Source the AmberTools environment first:

```bash
source ${AMBERHOME}/amber.sh
```

**`g16` not found.**
The Gaussian environment must be set up (`g16root`, profile script) so
that `g16` is on `PATH`. Make sure that `g16 < input.gjf > input.log`
works in your shell before running the script.

**RESP fails complaining about ESP data.**
If you passed a `gout` directly and RESP fails, the log may not actually
contain `ESP Fit Center` records. Re-run Gaussian with
`Pop=MK IOp(6/33=2,...)` or let the script generate and run the QM job.

**`parmchk2 -a Y` does not emit a `NONBON` section.**
Very old AmberTools releases may behave differently. The script raises a
clear error listing the missing atom types so that you can rebuild
AmberTools or fall back to a manual procedure.

**Output coordinates differ from input.**
Make sure your `g16` run used the script's generated input (which contains
`NoSymm`). A `gout` that was produced without `NoSymm` will have been
reoriented to Gaussian's standard orientation, and antechamber will read
those reoriented coordinates.

---

## License

MIT &mdash; see `LICENSE`.

## Citation

If you use this script in published work, please cite AmberTools and
Gaussian 16, as well as the GAFF / GAFF2 papers as appropriate.
