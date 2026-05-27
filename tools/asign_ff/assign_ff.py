#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
assign_ff.py
============

Interactive tool that assigns an empirical force field (GAFF / GAFF2, etc.) to
an arbitrary small molecule and writes out parameter / topology / structure
files in Amber (prep + frcmod), Gromacs (top + gro), and 3D-RISM (.rsm) formats.

External programs used (must be on PATH in the runtime environment):
    - antechamber, parmchk2   (AmberTools; always used)
    - tleap                   (AmberTools; only when Gromacs output is selected)
    - g16                     (Gaussian 16; only when RESP needs a QM calculation)
Python dependency:
    - parmed  (bundled with AmberTools; needed only for Gromacs output.
               Install with `pip install parmed` if missing.)

Input formats (auto-detected from extension + content):
    pdb, xyz, mol2, sdf (= mol), gjf (= com/gau, Gaussian input),
    gout (= log/out, Gaussian output)

Usage:
    python assign_ff.py molecule.mol2
    (If no input file is given as an argument, the path is asked for.)

Workflow overview:
    1) Enter a residue name (up to 3 chars, default MOL)
    2) Auto-detect the input format (xyz is not supported by antechamber,
       so it is converted to PDB internally)
    3) Enter total charge and spin multiplicity
    4) Select a force field (GAFF / GAFF2)
    5) Select a charge method (AM1-BCC / RESP)
    6) Charge calculation
         AM1-BCC : antechamber -c bcc
         RESP    : if the input is a gout that already contains ESP data,
                   read it directly with -c resp. Otherwise run g16 (asking
                   whether to optimize, default = no) and read RESP from its
                   output.
    7) Output (multiple selections allowed)
         Amber   : antechamber (prepi) + parmchk2 (frcmod, missing only)
         Gromacs : tleap builds prmtop/inpcrd -> ParmEd converts to top/gro
         RISM    : parmchk2 -a Y provides the full LJ table, joined with the
                   mol2 to write the .rsm (does not depend on tleap / ParmEd)
"""

import os
import sys
import re
import shlex
import shutil
import subprocess

# ===========================================================================
# Settings (edit as needed)
# ===========================================================================
# Theory level etc. for the RESP QM job. The GAFF/GAFF2 standard is HF/6-31G* + Pop=MK.
# NoSymm prevents Gaussian from reorienting to its "standard orientation", so the
# geometry antechamber reads back from the output keeps the original input frame.
QM_METHOD      = "HF/6-31G*"
QM_EXTRA       = "SCF=Tight NoSymm Pop=MK IOp(6/33=2,6/41=10,6/42=17)"
DEFAULT_DO_OPT = False         # Default for "optimize geometry in the RESP QM job?"
G16_NPROC      = 4             # Gaussian parallelism (%nprocshared)
G16_MEM        = "2GB"         # Gaussian memory (%mem)
G16_EXE        = "g16"         # Gaussian executable name

KCAL2J         = 4184.0        # 1 kcal = 4184 J (epsilon: kcal/mol -> J/mol)

RESNAME        = "MOL"         # Residue name; overwritten interactively at runtime.

# ===========================================================================
# Low-level utilities
# ===========================================================================
def die(msg):
    print(f"\n[ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


def need_exe(name):
    """Check that an external command exists. Exit if not. Return its path."""
    path = shutil.which(name)
    if path is None:
        die(f"Executable '{name}' not found. Please check your PATH.")
    return path


def run(cmd, cwd, stdin_path=None, stdout_path=None):
    """Run a subprocess (cmd is a list). On failure, print output and exit.

    stdin_path and stdout_path are resolved RELATIVE TO cwd, because the
    subprocess runs in cwd while Python's open() would otherwise use the
    interpreter's own working directory.
    """
    printable = " ".join(shlex.quote(c) for c in cmd)
    if stdin_path:
        printable += f" < {stdin_path}"
    if stdout_path:
        printable += f" > {stdout_path}"
    print(f"    $ {printable}")

    fin = open(os.path.join(cwd, stdin_path), "r") if stdin_path else None
    fout = open(os.path.join(cwd, stdout_path), "w") if stdout_path else None
    try:
        proc = subprocess.run(
            cmd, cwd=cwd, stdin=fin,
            stdout=(fout if fout else subprocess.PIPE),
            stderr=subprocess.PIPE, text=True,
        )
    finally:
        if fin:
            fin.close()
        if fout:
            fout.close()

    if proc.returncode != 0:
        if proc.stdout:
            print(proc.stdout, file=sys.stderr)
        if proc.stderr:
            print(proc.stderr, file=sys.stderr)
        die(f"Command exited abnormally (returncode={proc.returncode}):\n      {printable}")
    return proc


# ===========================================================================
# Input format auto-detection
# ===========================================================================
# Internal name -> antechamber -fi value
FI_OF = {
    "pdb":  "pdb",
    "mol2": "mol2",
    "sdf":  "mdl",     # In antechamber an MDL molfile = "mdl" (reads .sdf/.mol)
    "gjf":  "gcrt",    # Gaussian Cartesian input
    "gout": "gout",    # Gaussian output
    # xyz is not supported by antechamber -> convert to PDB, then treat as "pdb"
}

EXT_MAP = {
    ".pdb": "pdb", ".ent": "pdb",
    ".xyz": "xyz",
    ".mol2": "mol2",
    ".sdf": "sdf", ".mol": "sdf",
    ".gjf": "gjf", ".com": "gjf", ".gau": "gjf",
    ".gout": "gout", ".log": "gout", ".out": "gout",
}


def detect_format(path):
    """Prefer the extension; if ambiguous, inspect the contents."""
    ext = os.path.splitext(path)[1].lower()
    fmt = EXT_MAP.get(ext)

    try:
        with open(path, "r", errors="ignore") as f:
            head = f.read(8192)
    except Exception as e:
        die(f"Cannot read the input file: {e}")

    # Content-based detection / override
    if "@<TRIPOS>MOLECULE" in head:
        return "mol2"
    if ("Entering Gaussian System" in head or "Gaussian(R)" in head
            or "Normal termination of Gaussian" in head):
        return "gout"
    if re.search(r"V2000|V3000", head) and "M  END" in head:
        return "sdf"
    if "$$$$" in head:
        return "sdf"
    if re.search(r"(?m)^\s*(ATOM|HETATM|CRYST1)", head):
        return "pdb"
    # Gaussian input: has a %link0 line or a route line starting with '#'
    # (gout was already excluded above)
    if re.search(r"(?m)^\s*%(chk|mem|nproc)", head, re.I) or re.search(r"(?m)^\s*#", head):
        return "gjf"
    # xyz: line 1 is an integer atom count, line 2 a comment, then "element x y z"
    lines = head.splitlines()
    if len(lines) >= 3 and lines[0].strip().isdigit():
        toks = lines[2].split()
        if len(toks) >= 4:
            try:
                float(toks[1]); float(toks[2]); float(toks[3])
                return "xyz"
            except ValueError:
                pass

    if fmt:
        return fmt
    die("Could not determine the input format. Supported: pdb, xyz, mol2, sdf, gjf, gout")


# ===========================================================================
# xyz -> PDB conversion (antechamber cannot read xyz; bonds are perceived from
# interatomic distances)
# ===========================================================================
def xyz_to_pdb(xyz_path, pdb_path):
    with open(xyz_path) as f:
        lines = f.read().splitlines()
    try:
        natom = int(lines[0].split()[0])
    except (IndexError, ValueError):
        die("Could not parse the atom count on line 1 of the xyz file.")
    atom_lines = [ln for ln in lines[2:2 + natom] if ln.strip()]
    if len(atom_lines) < natom:
        die("Atom count and number of coordinate lines in the xyz file do not match.")

    counts = {}
    with open(pdb_path, "w") as out:
        for i, ln in enumerate(atom_lines, start=1):
            t = ln.split()
            elem = t[0].strip()
            elem = elem[0].upper() + elem[1:].lower() if len(elem) > 1 else elem.upper()
            x, y, z = float(t[1]), float(t[2]), float(t[3])
            counts[elem] = counts.get(elem, 0) + 1
            name = f"{elem}{counts[elem]}"[:4]
            # Standard PDB column alignment; element symbol also placed in cols 77-78.
            out.write(
                "HETATM{:>5d} {:<4s}{:1s}{:>3s} {:1s}{:>4d}    "
                "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n".format(
                    i, name, "", RESNAME, "A", 1, x, y, z, 1.00, 0.00, elem.upper()
                )
            )
        out.write("END\n")
    print(f"    Converted xyz to PDB: {os.path.basename(pdb_path)} "
          f"(bonds will be perceived by antechamber from interatomic distances).")


# ===========================================================================
# Decide whether a gout already contains ESP data for RESP
# ===========================================================================
def gout_has_esp(path):
    """Return True if the ESP fit points written by iop(6/33=2) are present.

    The appearance of 'ESP Fit Center' is the most reliable indicator.
    """
    try:
        with open(path, "r", errors="ignore") as f:
            txt = f.read()
    except Exception:
        return False
    if "ESP Fit Center" in txt:
        return True
    if "Electrostatic Properties Using The SCF Density" in txt and "Fit " in txt:
        return True
    return False


# ===========================================================================
# Center-of-mass translation of a charged mol2 (in place)
# ===========================================================================
# Atomic masses (u) keyed by element symbol. Used to compute the mass center.
ATOMIC_MASS = {
    "H": 1.008,  "He": 4.0026, "Li": 6.94,   "Be": 9.0122, "B": 10.81,
    "C": 12.011, "N": 14.007,  "O": 15.999,  "F": 18.998,  "Ne": 20.180,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06,  "Cl": 35.45,  "Ar": 39.948, "K": 39.098,  "Ca": 40.078,
    "Br": 79.904, "I": 126.90,
}


def element_from_mol2(atom_name, sybyl_type):
    """Best-effort element symbol from a TRIPOS atom (SYBYL type or atom name)."""
    base = sybyl_type.split(".")[0]
    if base in ATOMIC_MASS:
        return base
    cand = base[:1].upper() + base[1:2].lower()
    if cand in ATOMIC_MASS:
        return cand
    m = re.match(r"([A-Za-z]{1,2})", atom_name)
    if m:
        s = m.group(1)
        for c in (s[:1].upper() + s[1:2].lower(), s[:1].upper()):
            if c in ATOMIC_MASS:
                return c
    return base[:1].upper()


def translate_mol2_to_com(work, mol2):
    """Translate all atoms in a TRIPOS mol2 so the mass center sits at the origin.

    This is a pure translation (no rotation); only the X/Y/Z columns of the
    ATOM block are rewritten. Charges and connectivity are untouched.
    """
    path = os.path.join(work, mol2)
    with open(path, errors="ignore") as f:
        lines = f.readlines()

    # Locate the ATOM block.
    start = end = None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s == "@<TRIPOS>ATOM":
            start = i + 1
        elif start is not None and s.startswith("@<TRIPOS>"):
            end = i
            break
    if start is None:
        die("Could not find an ATOM block in the mol2 file for COM translation.")
    if end is None:
        end = len(lines)

    # First pass: accumulate mass-weighted coordinates.
    msum = mx = my = mz = 0.0
    for i in range(start, end):
        p = lines[i].split()
        if len(p) < 6:
            continue
        elem = element_from_mol2(p[1], p[5])
        mass = ATOMIC_MASS.get(elem, 12.011)
        x, y, z = float(p[2]), float(p[3]), float(p[4])
        msum += mass
        mx += mass * x
        my += mass * y
        mz += mass * z
    if msum == 0.0:
        die("Total mass evaluated to zero during COM translation.")
    cx, cy, cz = mx / msum, my / msum, mz / msum

    # Second pass: rewrite coordinates, preserving the rest of each line.
    # TRIPOS ATOM record: id name x y z type subst_id subst_name charge
    for i in range(start, end):
        raw = lines[i].rstrip("\n")
        p = raw.split()
        if len(p) < 6:
            continue
        x = float(p[2]) - cx
        y = float(p[3]) - cy
        z = float(p[4]) - cz
        rest = p[5:]
        lines[i] = (
            "{:>7s} {:<8s} {:>10.4f} {:>10.4f} {:>10.4f} {}\n".format(
                p[0], p[1], x, y, z, " ".join(rest)
            )
        )

    with open(path, "w") as f:
        f.writelines(lines)
    print(f"    Translated geometry so the mass center is at the origin "
          f"(shift = {cx:.4f}, {cy:.4f}, {cz:.4f}).")


# ===========================================================================
# Interactive prompts
# ===========================================================================
def ask_int(prompt):
    while True:
        s = input(prompt).strip()
        try:
            return int(s)
        except ValueError:
            print("    Please enter an integer.")


def ask_yesno(prompt, default=False):
    d = "Y/n" if default else "y/N"
    while True:
        s = input(f"{prompt} [{d}]: ").strip().lower()
        if s == "":
            return default
        if s in ("y", "yes"):
            return True
        if s in ("n", "no"):
            return False
        print("    Please answer y or n.")


def ask_resname():
    while True:
        s = input("Enter residue name (alphanumeric, up to 3 chars, default MOL): ").strip().upper()
        if s == "":
            return "MOL"
        if s.isalnum() and 1 <= len(s) <= 3:
            return s
        print("    Please enter up to 3 alphanumeric characters.")


def ask_choice(prompt, options):
    """options: [(key, label), ...] single selection. Returns the chosen key."""
    print(prompt)
    for i, (_, label) in enumerate(options, start=1):
        print(f"    [{i}] {label}")
    while True:
        s = input("    Enter a number: ").strip()
        if s.isdigit() and 1 <= int(s) <= len(options):
            return options[int(s) - 1][0]
        print("    Please enter a valid number.")


def ask_multi(prompt, options):
    """Multiple selection (comma separated). Returns the list of chosen keys."""
    print(prompt)
    for i, (_, label) in enumerate(options, start=1):
        print(f"    [{i}] {label}")
    while True:
        s = input("    Enter number(s) separated by commas (e.g. 1,3): ").strip()
        toks = [t.strip() for t in s.split(",") if t.strip()]
        if toks and all(t.isdigit() and 1 <= int(t) <= len(options) for t in toks):
            keys = []
            for t in toks:
                k = options[int(t) - 1][0]
                if k not in keys:
                    keys.append(k)
            return keys
        print("    Please enter valid number(s).")


# ===========================================================================
# Charge calculation with antechamber (produces a charged mol2)
# ===========================================================================
def run_charges(work, stem, in_name, in_fmt, ff_at, charge_method, nc, mult):
    """Create a charged mol2 ('<stem>.mol2') in work and return its name."""
    mol2 = f"{stem}.mol2"
    fi = FI_OF[in_fmt]

    if charge_method == "bcc":
        # ---- AM1-BCC: no Gaussian needed --------------------------------
        run([ANTECHAMBER, "-i", in_name, "-fi", fi,
             "-o", mol2, "-fo", "mol2",
             "-c", "bcc", "-nc", str(nc), "-m", str(mult),
             "-at", ff_at, "-rn", RESNAME, "-pf", "y", "-s", "2"], cwd=work)
        return mol2

    # ---- RESP -----------------------------------------------------------
    direct_gout = None
    if in_fmt == "gout" and gout_has_esp(os.path.join(work, in_name)):
        print("    Input gout already contains ESP data; reading RESP directly "
              "without running g16.")
        direct_gout = in_name

    if direct_gout is None:
        # g16 is required -> ask whether to optimize (default = no)
        need_exe(G16_EXE)
        do_opt = ask_yesno("    Run geometry optimization in the RESP QM calculation?",
                           default=DEFAULT_DO_OPT)
        route = "#P {} {}{}".format(QM_METHOD, "Opt " if do_opt else "", QM_EXTRA)
        gjf = f"{stem}_resp.gjf"
        gout = f"{stem}_resp.gout"
        run([ANTECHAMBER, "-i", in_name, "-fi", fi,
             "-o", gjf, "-fo", "gcrt",
             "-nc", str(nc), "-m", str(mult),
             "-gk", route,
             "-gn", f"%nprocshared={G16_NPROC}", "-gm", f"%mem={G16_MEM}",
             "-pf", "y", "-s", "2"], cwd=work)

        print("    Running Gaussian 16 QM calculation ({}ESP)...".format(
            "geometry optimization + " if do_opt else ""))
        run([G16_EXE], cwd=work, stdin_path=gjf, stdout_path=gout)
        with open(os.path.join(work, gout), errors="ignore") as f:
            if "Normal termination" not in f.read():
                die(f"Gaussian did not terminate normally. Check {gout}.")
        direct_gout = gout

    run([ANTECHAMBER, "-i", direct_gout, "-fi", "gout",
         "-o", mol2, "-fo", "mol2",
         "-c", "resp", "-nc", str(nc), "-m", str(mult),
         "-at", ff_at, "-rn", RESNAME, "-pf", "y", "-s", "2"], cwd=work)
    return mol2


# ===========================================================================
# frcmod / prep generation
# ===========================================================================
def run_parmchk_missing(work, stem, mol2, ff_name):
    """frcmod with only the missing parameters (for Amber output / tleap)."""
    frcmod = f"{stem}.frcmod"
    run([PARMCHK2, "-i", mol2, "-f", "mol2", "-o", frcmod, "-s", ff_name], cwd=work)
    return frcmod


def run_parmchk_all(work, stem, mol2, ff_name):
    """frcmod with ALL parameters including pre-existing ones (-a Y),
    used to obtain the full LJ table for RISM."""
    allf = f"{stem}_all.frcmod"
    run([PARMCHK2, "-i", mol2, "-f", "mol2", "-o", allf,
         "-s", ff_name, "-a", "Y"], cwd=work)
    return allf


def write_prep(work, stem, mol2):
    prep = f"{stem}.prep"
    run([ANTECHAMBER, "-i", mol2, "-fi", "mol2",
         "-o", prep, "-fo", "prepi", "-rn", RESNAME, "-pf", "y", "-s", "2"], cwd=work)
    return prep


# ===========================================================================
# Gromacs output (build prmtop with tleap, convert with ParmEd)
# ===========================================================================
def build_prmtop(work, stem, mol2, frcmod, ff_name):
    leaprc = "leaprc.gaff2" if ff_name == "gaff2" else "leaprc.gaff"
    prmtop, inpcrd = f"{stem}.prmtop", f"{stem}.inpcrd"
    leapin = f"{stem}_tleap.in"
    with open(os.path.join(work, leapin), "w") as f:
        f.write(f"source {leaprc}\n")
        f.write(f"{RESNAME} = loadmol2 {mol2}\n")
        f.write(f"loadamberparams {frcmod}\n")
        f.write(f"saveamberparm {RESNAME} {prmtop} {inpcrd}\n")
        f.write("quit\n")
    run([TLEAP, "-s", "-f", leapin], cwd=work)
    if not os.path.exists(os.path.join(work, prmtop)):
        die("tleap failed to generate the prmtop. Check leap.log.")
    return prmtop, inpcrd


def write_gromacs(work, stem, prmtop, inpcrd):
    try:
        import parmed as pmd
    except ImportError:
        die("Gromacs output requires ParmEd (pip install parmed).")
    top = f"{stem}_GMX.top"
    gro = f"{stem}_GMX.gro"
    struct = pmd.load_file(os.path.join(work, prmtop), xyz=os.path.join(work, inpcrd))
    struct.save(os.path.join(work, top), format="gromacs", overwrite=True)
    struct.save(os.path.join(work, gro), overwrite=True)
    return [top, gro]


# ===========================================================================
# RISM output (built from parmchk2 -a Y + mol2; no tleap / ParmEd dependency)
# ===========================================================================
def parse_mol2_atoms(path):
    """Read (name, type, charge, x, y, z) from the ATOM block of a TRIPOS mol2."""
    atoms = []
    in_atom = False
    with open(path, errors="ignore") as f:
        for ln in f:
            s = ln.strip()
            if s.startswith("@<TRIPOS>"):
                in_atom = (s == "@<TRIPOS>ATOM")
                continue
            if in_atom and s:
                p = s.split()
                name = p[1]
                x, y, z = float(p[2]), float(p[3]), float(p[4])
                atype = p[5]
                try:
                    charge = float(p[8])
                except (IndexError, ValueError):
                    charge = float(p[-1])
                atoms.append((name, atype, charge, x, y, z))
    if not atoms:
        die("Could not read atom information from the mol2 file.")
    return atoms


def parse_nonbon(frcmod_path):
    """Parse the NONBON section of an frcmod:
    type -> (rstar = Rmin/2 [Angstrom], epsilon [kcal/mol])."""
    lj = {}
    in_nb = False
    section_kw = {"MASS", "BOND", "ANGLE", "ANGL", "DIHE", "IMPROPER",
                  "IMPROP", "NONBON", "CMAP", "LJEDIT"}
    with open(frcmod_path, errors="ignore") as f:
        for ln in f:
            tok = ln.split()
            if not tok:
                if in_nb:
                    break                 # NONBON ends at a blank line
                continue
            if tok[0].upper() in section_kw:
                in_nb = (tok[0].upper() == "NONBON")
                continue
            if in_nb and len(tok) >= 3:
                try:
                    lj[tok[0]] = (float(tok[1]), float(tok[2]))
                except ValueError:
                    pass
    return lj


def write_rism(work, stem, mol2, ff_name, label):
    """Write the 3D-RISM .rsm file according to the spec:

        line 1     : " $UDATA"                 (one leading space)
        line 2     : natom (10-digit integer) label
        line 3..   : atom_label  sigma[Angstrom]  epsilon[J/mol]  charge[e]  X  Y  Z
        last line  : " $END"                    (one leading space)

    sigma and epsilon are taken from the full NONBON table written by
    parmchk2 -a Y (including pre-existing parameters); charge and coordinates
    are taken from the charged mol2. This path does not depend on tleap.
    """
    allf = run_parmchk_all(work, stem, mol2, ff_name)
    atoms = parse_mol2_atoms(os.path.join(work, mol2))
    lj = parse_nonbon(os.path.join(work, allf))

    missing = sorted({a[1] for a in atoms if a[1] not in lj})
    if missing:
        die("Could not obtain NONBON parameters for atom type(s): "
            + ", ".join(missing)
            + f"\n      Check the NONBON section of {allf}.")

    rsm = f"{stem}.rsm"
    with open(os.path.join(work, rsm), "w") as f:
        f.write(" $UDATA\n")
        f.write(f"{len(atoms):10d} {label}\n")
        for (name, atype, charge, x, y, z) in atoms:
            rstar, eps_kcal = lj[atype]
            sigma = rstar * 2.0 ** (5.0 / 6.0)    # Rmin/2 -> sigma [Angstrom]
            eps_jmol = eps_kcal * KCAL2J           # kcal/mol -> J/mol
            f.write(
                "{:<6s} {:14.6f} {:16.6f} {:14.8f} "
                "{:14.6f} {:14.6f} {:14.6f}\n".format(
                    name, sigma, eps_jmol, charge, x, y, z
                )
            )
        f.write(" $END\n")
    return [rsm]


# ===========================================================================
# Main
# ===========================================================================
def main():
    global ANTECHAMBER, PARMCHK2, TLEAP, RESNAME

    # --- Input file ---
    if len(sys.argv) >= 2:
        in_path = sys.argv[1]
    else:
        in_path = input("Enter path to the input file: ").strip()
    in_path = os.path.abspath(in_path)
    if not os.path.isfile(in_path):
        die(f"Input file does not exist: {in_path}")

    stem = os.path.splitext(os.path.basename(in_path))[0]   # without extension = RISM label
    in_fmt = detect_format(in_path)
    print(f"\nInput file     : {in_path}")
    print(f"Detected format: {in_fmt}")

    # --- Always-required external commands ---
    ANTECHAMBER = need_exe("antechamber")
    PARMCHK2    = need_exe("parmchk2")
    TLEAP       = None   # resolved only when Gromacs output is selected

    # --- Residue name (asked first) ---
    print()
    RESNAME = ask_resname()

    # --- Prepare working directory ---
    work = os.path.abspath(f"{stem}_param")
    os.makedirs(work, exist_ok=True)

    if in_fmt == "xyz":
        in_name = f"{stem}_from_xyz.pdb"
        xyz_to_pdb(in_path, os.path.join(work, in_name))
        in_fmt = "pdb"
    else:
        in_name = os.path.basename(in_path)
        shutil.copy(in_path, os.path.join(work, in_name))

    # --- Total charge and spin multiplicity ---
    print()
    nc   = ask_int("Enter total molecular charge (integer): ")
    mult = ask_int("Enter spin multiplicity (integer, e.g. 1 = singlet, 2 = doublet): ")

    # --- Force field ---
    print()
    ff = ask_choice("Select the force field to use:",
                    [("gaff",  "GAFF  (General Amber Force Field)"),
                     ("gaff2", "GAFF2 (General Amber Force Field v2)")])

    # --- Charge method ---
    print()
    cm = ask_choice("Select the charge method:",
                    [("bcc",  "AM1-BCC (fast, no QM required)"),
                     ("resp", "RESP    (HF/6-31G* ESP, uses Gaussian 16)")])

    # --- Output format(s) (multiple allowed) ---
    print()
    outs = ask_multi("Select output format(s) (multiple allowed):",
                     [("amber",   "Amber   (prep + frcmod)"),
                      ("gromacs", "Gromacs (top + gro)"),
                      ("rism",    "3D-RISM (.rsm)")])

    # --- Geometry frame for the output coordinates ---
    print()
    geom = ask_choice("Select the geometry for the output coordinates:",
                      [("orig", "Keep the original input geometry"),
                       ("com",  "Translate the mass center to the origin")])

    # --- Charge calculation -> mol2 ---
    print("\n[Charge calculation and atom typing (antechamber)]")
    mol2 = run_charges(work, stem, in_name, in_fmt, ff, cm, nc, mult)

    # --- Optional: move mass center to the origin (pure translation) ---
    # NoSymm in the QM route already prevents Gaussian reorientation, so by
    # default the original input frame is preserved across all output formats.
    if geom == "com":
        translate_mol2_to_com(work, mol2)

    produced = []

    # --- Amber: prep + frcmod ---
    if "amber" in outs:
        print("\n[Amber output (prep + frcmod)]")
        frcmod = run_parmchk_missing(work, stem, mol2, ff)
        prep = write_prep(work, stem, mol2)
        produced += [prep, frcmod]

    # --- Gromacs: tleap -> ParmEd ---
    if "gromacs" in outs:
        print("\n[Gromacs output (top + gro)]")
        TLEAP = need_exe("tleap")
        gmx_frcmod = run_parmchk_missing(work, stem, mol2, ff)
        prmtop, inpcrd = build_prmtop(work, stem, mol2, gmx_frcmod, ff)
        produced += write_gromacs(work, stem, prmtop, inpcrd)

    # --- RISM: parmchk2 -a Y + mol2 ---
    if "rism" in outs:
        print("\n[3D-RISM output (.rsm)]")
        produced += write_rism(work, stem, mol2, ff, label=stem)

    # --- Summary ---
    print("\n=== Done ===")
    print(f"Residue name    : {RESNAME}")
    print(f"Output directory: {work}")
    for fn in produced:
        p = os.path.join(work, fn)
        mark = "OK " if os.path.exists(p) else "?? "
        print(f"  [{mark}] {fn}")


if __name__ == "__main__":
    main()
