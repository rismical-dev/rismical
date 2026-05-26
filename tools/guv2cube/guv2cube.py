#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
guv2cube.py
===========
Convert a 3D-RISM .guv distribution-function file into a Gaussian .cube file.
If a structure file (.xyz / .pdb / .mol2) is supplied, its atoms are embedded
into the cube.

Usage
-----
    python3 guv2cube.py <input.guv> [structure] [out.cub]

  Arguments after the first are classified automatically by extension
  (order-independent):
    .xyz / .pdb / .mol2 / .ent  -> structure file to embed
    .cub / .cube                -> output file name
    anything else               -> treated as the output file name
  If the output name is omitted, it defaults to <input>.cub.
  Examples:
    python3 guv2cube.py water.guv                 # no atoms (dummy atom if multi-map)
    python3 guv2cube.py water.guv solute.xyz      # embed atoms from solute.xyz
    python3 guv2cube.py water.guv solute.pdb g.cub

Input .guv format
-----------------
Output of the Fortran subroutine `write3dfunc`; ASCII / BIN is auto-detected.

  ASCII : line 1 = "##  " + comment, line 2 (format 4i8,6f16.8) =
          nvuq, N, N, N, rn, rn, rn, shift, shift, shift (rn = N*rdelta3d = box edge [Ang]),
          followed by nvuq*N**3 values (format e16.8e3). Order: iv outer, ig inner.
  BIN   : direct access (recl in bytes). bytes 0-15 = int32 x4 (nvuq,N,N,N),
          bytes 16-39 = rn x3, bytes 40-63 = shift x3, bytes 64- = float64 data.
          Endianness is auto-detected.

func3d is dimensionless; rdelta3d is in Angstrom. Cube coordinates are written
in Bohr (Gaussian default).

Conventions / assumptions (configurable via the flags below)
------------------------------------------------------------
* Grid linearization defaults to "x fastest": ig = ix + iy*N + iz*N^2
  (set SOURCE_FASTEST="z" for z fastest).
* The grid is assumed origin-centered: the corner (first grid point) is placed
  at shift - box/2 (set CENTER_GRID=False to place the corner at shift).
* Structure coordinates are assumed to share the RISM grid frame and are
  embedded as-is (Ang -> Bohr). Use STRUCT_SHIFT_ANG to translate if needed.
* nvuq==1 -> single-dataset cube. nvuq>1 -> multiple datasets in one file
  (requires a negative atom count; uses the structure atoms if given, otherwise
  one dummy atom). For viewers that cannot read multi-dataset cubes (e.g. VESTA),
  set SPLIT=True to write one file per map.
"""

import sys
import struct
import numpy as np

# ---- Configurable convention flags --------------------------------------
SOURCE_FASTEST  = "x"           # "x": ig=ix+iy*N+iz*N^2 / "z": ig=iz+iy*N+ix*N^2
CENTER_GRID     = True          # True: corner = shift-box/2 / False: corner = shift
SPLIT           = False         # True: split nvuq>1 into one file per map
DUMMY_Z         = 1             # atomic number of the dummy atom for multi-dataset cubes with no structure
STRUCT_SHIFT_ANG = (0.0, 0.0, 0.0)  # translation [Ang] applied to the structure (for alignment)
# -------------------------------------------------------------------------

BOHR = 0.52917721090380   # 1 Bohr in Angstrom

# Element symbol -> atomic number
_ELEMENTS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co "
    "Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb "
    "Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os "
    "Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md "
    "No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()
SYMBOL2Z = {s.upper(): i + 1 for i, s in enumerate(_ELEMENTS)}


# =========================================================================
#  Reading the .guv file
# =========================================================================
def is_ascii(path):
    with open(path, "rb") as f:
        return f.read(2) == b"##"


def read_ascii(path):
    with open(path) as f:
        title = f.readline().rstrip("\n")
        hdr_line = f.readline()
        toks = hdr_line.split()
        if len(toks) >= 10:
            nvuq = int(toks[0]); N = int(toks[1])
            box  = (float(toks[4]), float(toks[5]), float(toks[6]))
            shift = (float(toks[7]), float(toks[8]), float(toks[9]))
        else:
            # Fallback for fixed-width fields that run together (4i8,6f16.8)
            s = hdr_line.rstrip("\n")
            nvuq = int(s[0:8]); N = int(s[8:16])
            box  = (float(s[32:48]), float(s[48:64]), float(s[64:80]))
            shift = (float(s[80:96]), float(s[96:112]), float(s[112:128]))
        ng3d = N ** 3
        vals = np.array(f.read().split(), dtype=np.float64)
    if vals.size != nvuq * ng3d:
        raise ValueError(f"value count mismatch: read {vals.size}, expected {nvuq*ng3d}")
    return nvuq, N, box, shift, vals.reshape(nvuq, ng3d), title


def read_bin(path):
    raw = open(path, "rb").read()
    for endian in ("<", ">"):
        if len(raw) < 64:
            break
        nvuq, n1, n2, n3 = struct.unpack(endian + "4i", raw[:16])
        if 0 < nvuq < 1_000_000 and 0 < n1 < 100_000 and n1 == n2 == n3:
            ng3d = n1 ** 3
            need = 64 + 8 * nvuq * ng3d
            if len(raw) >= need:
                box   = struct.unpack(endian + "3d", raw[16:40])
                shift = struct.unpack(endian + "3d", raw[40:64])
                data = np.frombuffer(raw, dtype=np.dtype(endian + "f8"),
                                     count=nvuq * ng3d, offset=64).reshape(nvuq, ng3d).copy()
                tag = "little-endian" if endian == "<" else "big-endian"
                return nvuq, n1, box, shift, data, f"RISM guv (binary, {tag})"
    raise ValueError("cannot parse binary .guv (endianness/size mismatch)")


# =========================================================================
#  Reading structure files (.xyz / .pdb / .mol2) -> [(Z, x, y, z)] in Angstrom
# =========================================================================
def _sym_to_z(token):
    """Element symbol or atomic-number string -> atomic number. None if unknown."""
    token = token.strip()
    if not token:
        return None
    if token.lstrip("+-").isdigit():          # written as an atomic number
        z = int(token)
        return z if 1 <= z <= 118 else None
    return SYMBOL2Z.get(token.capitalize().upper())


def _elem_from_name(name):
    """Infer the element symbol from an atom name (e.g. 'CA', 'OG1', 'Cl-')."""
    letters = "".join(c for c in name if c.isalpha())
    if not letters:
        return ""
    one = letters[0].upper()
    two = (letters[0].upper() + letters[1].lower()) if len(letters) >= 2 else ""
    # Prefer common single-letter elements (avoid reading 'CA' as Ca)
    if one in ("C", "N", "O", "H", "S", "P", "F", "B", "K", "I", "V", "W", "Y", "U"):
        return one
    if two in SYMBOL2Z:
        return two
    if one in SYMBOL2Z:
        return one
    return two


def _read_xyz(path):
    atoms = []
    lines = open(path).read().splitlines()
    try:
        n = int(lines[0].split()[0]); start = 2
    except (ValueError, IndexError):
        n = None; start = 0
    body = lines[start:start + n] if n is not None else lines[start:]
    for ln in body:
        t = ln.split()
        if len(t) < 4:
            continue
        z = _sym_to_z(t[0])
        if z is None:
            continue
        atoms.append((z, float(t[1]), float(t[2]), float(t[3])))
    return atoms


def _read_pdb(path):
    atoms = []
    for ln in open(path):
        if ln[:6] in ("ATOM  ", "HETATM"):
            try:
                x = float(ln[30:38]); y = float(ln[38:46]); zc = float(ln[46:54])
            except ValueError:
                continue
            elem = ln[76:78].strip()           # element symbol in columns 77-78 (standard)
            if not elem:
                elem = _elem_from_name(ln[12:16])
            z = _sym_to_z(elem)
            if z is None:
                continue
            atoms.append((z, x, y, zc))
    return atoms


def _read_mol2(path):
    atoms = []
    section = None
    for ln in open(path):
        s = ln.strip()
        if s.startswith("@<TRIPOS>"):
            section = s.upper()
            continue
        if section == "@<TRIPOS>ATOM" and s:
            t = s.split()
            if len(t) < 6:
                continue
            try:
                x, y, zc = float(t[2]), float(t[3]), float(t[4])
            except ValueError:
                continue
            z = _sym_to_z(t[5].split(".")[0])  # SYBYL atom type "C.3" -> "C"
            if z is None:
                z = _sym_to_z(_elem_from_name(t[1]))
            if z is None:
                continue
            atoms.append((z, x, y, zc))
    return atoms


def read_structure(path):
    ext = path.lower().rsplit(".", 1)[-1]
    if ext == "xyz":
        atoms = _read_xyz(path)
    elif ext in ("pdb", "ent"):
        atoms = _read_pdb(path)
    elif ext == "mol2":
        atoms = _read_mol2(path)
    else:
        raise ValueError(f"unsupported structure format: .{ext} (only xyz/pdb/mol2)")
    if not atoms:
        print(f"[guv2cube] warning: no atoms could be read from {path}.")
    return atoms


# =========================================================================
#  Grid reordering: ig -> cube output order [ix, iy, iz] (iz fastest)
# =========================================================================
def to_vol(flat1d, N):
    v = np.asarray(flat1d).reshape(N, N, N)   # C order -> last axis is fastest
    if SOURCE_FASTEST == "x":
        v = v.transpose(2, 1, 0)              # [iz,iy,ix] -> [ix,iy,iz]
    return np.ascontiguousarray(v)


# =========================================================================
#  Writing the cube file
# =========================================================================
def _write_data(o, flat):
    strs = np.char.mod("%13.5E", flat.astype(np.float64))
    for i in range(0, strs.size, 6):
        o.write("".join(strs[i:i + 6]) + "\n")


def write_one_cube(fname, vols, N, origin, d, title1, title2, atoms, multi):
    """vols: list of [ix,iy,iz] arrays. atoms: [(Z, charge, x, y, z)] in Bohr."""
    nat = len(atoms)
    natoms_field = -nat if multi else nat
    with open(fname, "w") as o:
        o.write(title1 + "\n")
        o.write(title2 + "\n")
        o.write(f"{natoms_field:5d}{origin[0]:12.6f}{origin[1]:12.6f}{origin[2]:12.6f}\n")
        o.write(f"{N:5d}{d[0]:12.6f}{0.0:12.6f}{0.0:12.6f}\n")
        o.write(f"{N:5d}{0.0:12.6f}{d[1]:12.6f}{0.0:12.6f}\n")
        o.write(f"{N:5d}{0.0:12.6f}{0.0:12.6f}{d[2]:12.6f}\n")
        for (z, ch, x, y, zc) in atoms:
            o.write(f"{z:5d}{ch:12.6f}{x:12.6f}{y:12.6f}{zc:12.6f}\n")
        if multi:
            m = len(vols)
            ids = [m] + list(range(1, m + 1))   # m, 1, 2, ..., m
            for i in range(0, len(ids), 10):     # format 10I5
                o.write("".join(f"{v:5d}" for v in ids[i:i + 10]) + "\n")
            flat = np.stack(vols, axis=-1).ravel(order="C")  # iv fastest, then iz,iy,ix
        else:
            flat = vols[0].ravel(order="C")                  # iz fastest
        _write_data(o, flat)


# =========================================================================
#  Main routine
# =========================================================================
def convert(inp, out=None, struct_file=None):
    if out is None:
        out = inp.rsplit(".", 1)[0] + ".cub"

    ascii_mode = is_ascii(inp)
    reader = read_ascii if ascii_mode else read_bin
    nvuq, N, box, shift, data, title = reader(inp)

    # Grid spacing [Ang] and corner coordinate [Ang] -> Bohr
    d_ang      = [box[k] / N for k in range(3)]
    corner_ang = [shift[k] - (box[k] / 2.0 if CENTER_GRID else 0.0) for k in range(3)]
    d_bohr      = [v / BOHR for v in d_ang]
    origin_bohr = [v / BOHR for v in corner_ang]
    center_bohr = [shift[k] / BOHR for k in range(3)]

    # Structure file -> cube atom rows (Bohr)
    real_atoms = []
    if struct_file:
        for (z, x, y, zc) in read_structure(struct_file):
            xb = (x + STRUCT_SHIFT_ANG[0]) / BOHR
            yb = (y + STRUCT_SHIFT_ANG[1]) / BOHR
            zb = (zc + STRUCT_SHIFT_ANG[2]) / BOHR
            real_atoms.append((z, float(z), xb, yb, zb))

    vols = [to_vol(data[iv], N) for iv in range(nvuq)]

    t1 = f"Gaussian cube from RISM guv  [{title.strip()}]"
    src = f"struct={struct_file}" if struct_file else "no-struct"
    t2 = (f"nvuq={nvuq} N={N} spacing={d_ang[0]:.6f}Ang box={box[0]:.6f}Ang "
          f"({'ASCII' if ascii_mode else 'BIN'}; {src})")

    written = []
    if nvuq == 1:
        write_one_cube(out, vols, N, origin_bohr, d_bohr, t1, t2,
                       atoms=real_atoms, multi=False)
        written.append(out)
    elif SPLIT:
        base = out[:-4] if out.lower().endswith(".cub") else out
        for iv in range(nvuq):
            f = f"{base}_{iv + 1}.cub"
            write_one_cube(f, [vols[iv]], N, origin_bohr, d_bohr,
                           t1 + f"  map {iv + 1}/{nvuq}", t2,
                           atoms=real_atoms, multi=False)
            written.append(f)
    else:
        # A multi-dataset cube requires a negative atom count
        atoms = real_atoms if real_atoms else \
            [(DUMMY_Z, 0.0, center_bohr[0], center_bohr[1], center_bohr[2])]
        write_one_cube(out, vols, N, origin_bohr, d_bohr, t1, t2,
                       atoms=atoms, multi=True)
        written.append(out)

    # Summary
    print(f"[guv2cube] input        : {inp}  ({'ASCII' if ascii_mode else 'BIN'})")
    print(f"[guv2cube] nvuq         : {nvuq}  (number of density maps)")
    print(f"[guv2cube] ngrid3d      : {N}  ({N**3} points/map)")
    print(f"[guv2cube] spacing[Ang] : {d_ang[0]:.6f}   box[Ang]: {box[0]:.6f}")
    print(f"[guv2cube] spacing[Bohr]: {d_bohr[0]:.6f}  origin[Bohr]: "
          f"({origin_bohr[0]:.4f}, {origin_bohr[1]:.4f}, {origin_bohr[2]:.4f})")
    if struct_file:
        print(f"[guv2cube] structure    : {struct_file}  ({len(real_atoms)} atoms)")
    for f in written:
        print(f"[guv2cube] output       : {f}")
    return written


def main(argv):
    if len(argv) < 2:
        print("usage: python3 guv2cube.py <input.guv> [structure.(xyz|pdb|mol2)] [out.cub]")
        # Default behavior (backward compatible)
        convert("sample.guv", "sample.cub")
        return
    inp = argv[1]
    out = None
    struct_file = None
    for a in argv[2:]:
        ext = a.lower().rsplit(".", 1)[-1] if "." in a else ""
        if ext in ("xyz", "pdb", "ent", "mol2"):
            struct_file = a
        else:                       # .cub/.cube or anything else: treat as output name
            out = a
    convert(inp, out, struct_file)


if __name__ == "__main__":
    main(sys.argv)
