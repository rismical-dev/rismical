# guv2cube

Convert 3D-RISM distribution-function files (`.guv`) into Gaussian **cube**
files (`.cub`), optionally embedding a solute structure
(`.xyz` / `.pdb` / `.mol2`) so the molecule and the distribution can be
visualized together.

The input `.guv` files are those written by the Fortran subroutine
`write3dfunc`. Both its **ASCII** and **binary** output formats are supported
and are detected automatically.

---

## Requirements

- Python 3.6+
- [NumPy](https://numpy.org/)

```bash
pip install numpy
```

---

## Usage

```bash
python3 guv2cube.py <input.guv> [structure] [out.cub]
```

The first argument is always the input `.guv` file. Any further arguments are
classified automatically **by file extension**, so the structure file and the
output name may be given in any order:

| Argument extension          | Interpreted as            |
|-----------------------------|---------------------------|
| `.xyz` `.pdb` `.mol2` `.ent`| structure file to embed   |
| `.cub` `.cube`              | output file name          |
| anything else               | output file name          |

If the output name is omitted, it defaults to `<input>.cub`.

### Examples

```bash
# Minimal: no structure embedded
python3 guv2cube.py water.guv                  # -> water.cub

# Embed a solute from an XYZ file
python3 guv2cube.py water.guv solute.xyz        # -> water.cub (with atoms)

# Embed a solute and choose the output name
python3 guv2cube.py water.guv solute.pdb g.cub  # -> g.cub

# Arguments are order-independent
python3 guv2cube.py water.guv g.cub solute.mol2 # -> g.cub
```

---

## Input `.guv` format

The script reads the output of `write3dfunc` and auto-detects the format by
inspecting the first two bytes (`##` => ASCII, otherwise binary).

**ASCII**
- Line 1: `"##  "` followed by an 80-character comment.
- Line 2 (Fortran format `4i8,6f16.8`):
  `nvuq, N, N, N, rn, rn, rn, shift, shift, shift`,
  where `rn = N * rdelta3d` is the box edge length in Å.
- Remaining lines: `nvuq * N**3` values (format `e16.8e3`),
  ordered with the map index `iv` outermost and the grid index `ig` innermost.

**Binary** (direct access, record length in bytes)
- bytes 0–15: four `int32` values `nvuq, N, N, N`
- bytes 16–39: three `float64` values `rn, rn, rn`
- bytes 40–63: three `float64` values `shift, shift, shift`
- bytes 64–: `nvuq * N**3` `float64` data values (same `iv`-outer ordering)
- Endianness (little / big) is detected automatically.

### Quantities

| Symbol     | Meaning                          | Unit        |
|------------|----------------------------------|-------------|
| `func3d`   | density / distribution value     | dimensionless |
| `nvuq`     | number of density maps           | –           |
| `ngrid3d`  | grid points per axis (`N`)       | –           |
| `rdelta3d` | grid spacing (`rn / N`)          | Å           |

Grid coordinates are converted from Å to **Bohr** on output, which is the
Gaussian cube default.

---

## Output cube format

The number of maps in the `.guv` file controls how the cube is written:

- **`nvuq == 1`** — a standard single-dataset cube
  (positive atom count in the header).
- **`nvuq > 1`** — a single multi-dataset cube. The Gaussian cube
  specification requires a **negative atom count** for multiple datasets, plus
  an extra "dataset IDs" line; the script writes both. Each grid point then
  stores all map values consecutively. GaussView and VMD can switch between the
  maps.

If a structure file is supplied, its atoms are embedded (atomic number, nuclear
charge, and Bohr coordinates). If no structure is given **and** `nvuq > 1`, a
single dummy atom is inserted at the box center, because the negative-atom-count
format requires at least one atom row.

### Structure file parsing

| Format  | Coordinates                                   | Element source                                   |
|---------|-----------------------------------------------|--------------------------------------------------|
| `.xyz`  | columns 2–4                                   | symbol (or atomic number) in column 1            |
| `.pdb`  | columns 31–38 / 39–46 / 47–54 of `ATOM`/`HETATM` | element field (cols 77–78); falls back to atom name |
| `.mol2` | `@<TRIPOS>ATOM` section                        | SYBYL atom type (e.g. `C.3` → C); falls back to atom name |

All structure coordinates are read in Å and converted to Bohr for the cube.

---

## Configuration flags

These constants near the top of `guv2cube.py` change conventions that cannot be
inferred from the file alone:

| Flag               | Default         | Effect                                                                 |
|--------------------|-----------------|------------------------------------------------------------------------|
| `SOURCE_FASTEST`   | `"x"`           | Grid linearization. `"x"`: `ig = ix + iy*N + iz*N**2`; `"z"`: z fastest. |
| `CENTER_GRID`      | `True`          | `True`: grid corner at `shift - box/2` (origin-centered); `False`: corner at `shift`. |
| `SPLIT`            | `False`         | If `True`, write one file per map (`name_1.cub`, `name_2.cub`, …) for viewers that cannot read multi-dataset cubes (e.g. VESTA). |
| `DUMMY_Z`          | `1`             | Atomic number of the dummy atom used for multi-dataset cubes when no structure is given. |
| `STRUCT_SHIFT_ANG` | `(0.0,0.0,0.0)` | Translation (Å) applied to the structure for alignment with the grid.  |

### Notes on the assumptions

- If the visualized distribution appears mirrored or rotated relative to the
  solute, switch `SOURCE_FASTEST` between `"x"` and `"z"`. For symmetric
  distributions on a cubic grid the two give identical results.
- The structure is assumed to share the RISM grid's coordinate frame and is
  embedded as-is. If the molecule and the density are offset, use
  `STRUCT_SHIFT_ANG` to translate it.

---

## Viewing the result

- **GaussView / VMD** — read single- and multi-dataset cubes; for `nvuq > 1`
  you can toggle between maps.
- **VESTA** — does not read multi-dataset cubes. Set `SPLIT = True` to produce
  one standard cube per map; the structure is embedded in each.

---

## Troubleshooting

- *"値の数が不一致 / value count mismatch"* — the data count does not match
  `nvuq * N**3`; the file is likely truncated or the wrong format.
- *"バイナリ .guv を解釈できません / cannot parse binary .guv"* — header values
  or file size are inconsistent with the expected layout (check that the file is
  a genuine `write3dfunc` binary output).
- *"未対応の構造形式 / unsupported structure format"* — only `.xyz`, `.pdb`,
  and `.mol2` (and `.ent`) are recognized as structures.
- Atoms missing from the cube — the parser could not assign an element; check
  the element/atom-name fields in the structure file.
