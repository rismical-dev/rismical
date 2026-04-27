import sys
import os
from pathlib import Path


def find_include_file(name, base_dir):
    """Search for include files"""
    # 1. Current directory (same level as the calling top file)
    p = Path(base_dir) / name
    if p.exists():
        return p

    # 2. GMXDATA/top (environment variable)
    gmx = os.environ.get("GMXDATA")
    if gmx:
        p = Path(gmx) / "top" / name
        if p.exists():
            return p

    # 3. Default path for Homebrew (fallback)
    p = Path("/opt/homebrew/share/gromacs/top") / name
    if p.exists():
        return p

    raise FileNotFoundError(f"include file not found: {name}")


def read_atomtypes_recursive(top_file, visited=None):
    """
    Read atomtypes from a top file and all included files recursively, returning a dictionary of atomtype to (sigma, epsilon).
    """
    if visited is None:
        visited = set()

    atom_data = {}
    top_file = Path(top_file).resolve()

    if top_file in visited:
        return atom_data
    visited.add(top_file)

    base_dir = top_file.parent
    atoms_section = False

    with open(top_file) as f:
        for line in f:
            line = line.strip()

            # ---- include processing ----
            if line.startswith("#include"):
                parts = line.split('"')
                if len(parts) >= 3:
                    name = parts[1]
                    try:
                        inc_file = find_include_file(name, base_dir)
                        atom_data.update(
                            read_atomtypes_recursive(inc_file, visited)
                        )
                    except FileNotFoundError:
                        # Skip missing include files with a warning, but continue processing other files
                        print(f"[Warning] Skipped missing include file: {name}", file=sys.stderr)
                continue

            # ---- atomtypes ----
            if line.startswith("[ atomtypes ]"):
                atoms_section = True
                continue

            if atoms_section:
                if line.startswith("["):
                    atoms_section = False
                    continue

                # Exclude inline comments (everything after ;) and skip empty lines
                line_no_comment = line.split(';')[0].strip()
                if not line_no_comment:
                    continue

                parts = line_no_comment.split()
                if len(parts) >= 6:
                    atomtype = parts[0]
                    try:
                        # Force field dependent column count, so take the last two columns (sigma, epsilon)
                        sigma = float(parts[-2]) * 10.0      # nm -> Angstrom
                        epsilon = float(parts[-1]) * 1000.0  # kJ/mol -> J/mol
                        atom_data[atomtype] = (sigma, epsilon)
                    except ValueError:
                        continue

    return atom_data


def read_gro(gro_file):
    coords = []

    with open(gro_file) as f:
        lines = f.readlines()

    n_atoms = int(lines[1].strip())

    for i in range(2, 2 + n_atoms):

        line = lines[i]

        x = float(line[20:28]) * 10.0 # nm -> Angstrom
        y = float(line[28:36]) * 10.0
        z = float(line[36:44]) * 10.0

        coords.append((x, y, z))

    return n_atoms, coords


def read_atoms_recursive(top_file, visited=None):
    if visited is None:
        visited = set()

    atoms = []
    top_file = Path(top_file).resolve()

    if top_file in visited:
        return atoms
    visited.add(top_file)

    base_dir = top_file.parent
    atoms_section = False

    with open(top_file) as f:
        for line in f:
            line = line.strip()

            # ---- include ----
            if line.startswith("#include"):
                parts = line.split('"')
                if len(parts) >= 3:
                    name = parts[1]
                    try:
                        inc_file = find_include_file(name, base_dir)
                        atoms.extend(
                            read_atoms_recursive(inc_file, visited)
                        )
                    except FileNotFoundError:
                        # Warning is already output in read_atomtypes_recursive, so just ignore here
                        pass
                continue

            # ---- atoms ----
            if line.startswith("[ atoms ]"):
                atoms_section = True
                continue

            if atoms_section:
                if line.startswith("["):
                    atoms_section = False
                    continue

                # Exclude inline comments (everything after ;) and skip empty lines
                line_no_comment = line.split(';')[0].strip()
                if not line_no_comment:
                    continue

                if line_no_comment[0].isdigit():
                    parts = line_no_comment.split()
                    nr = int(parts[0])
                    type_name = parts[1]
                    # Obtain atom_label as the 5th column (index 4) of [ atoms ].
                    # If the column is missing, use type_name as a fallback.
                    atom_label = parts[4] if len(parts) > 4 else type_name
                    charge = float(parts[6])
                    atoms.append((nr, type_name, atom_label, charge))

    return atoms


def main(top_file, gro_file):

    atomtypes = read_atomtypes_recursive(top_file)
    atoms = read_atoms_recursive(top_file)
    n_atoms, coords = read_gro(gro_file)

    # 1st line: number_of_atoms name_of_molecule(optional)
    molecule_name = Path(top_file).stem
    print(f"{n_atoms} {molecule_name}")

    # 2nd line and beyond: atom_label LJ_sigma[Angs] LJ_epsilon[J/mol] point_charge[e] coordinate_x, y, z [Angs]
    for i, (nr, type_name, atom_label, charge) in enumerate(atoms[:n_atoms]):

        if type_name not in atomtypes:
            raise KeyError(f"atomtype '{type_name}' not found in any include")

        sigma, epsilon = atomtypes[type_name]
        x, y, z = coords[i]

        # atom_label LJ_sigma[Angs] LJ_epsilon[J/mol] point_charge[e] coordinate_x, y, z [Angs]
        print(f"{atom_label:<8s} {sigma:8.4f} {epsilon:10.4f} {charge:8.4f} {x:8.4f} {y:8.4f} {z:8.4f}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py file.top file.gro")
    else:
        main(sys.argv[1], sys.argv[2])