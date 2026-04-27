import sys
import os
import argparse
from pathlib import Path

def find_include_file(name, base_dir):
    """Search for include files"""
    # 1. Current directory (same level as the calling top file)
    p = Path(base_dir) / name
    if p.exists():
        return p

    # 2. GMXDATA/top (Environment variable)
    gmx = os.environ.get("GMXDATA")
    if gmx:
        p = Path(gmx) / "top" / name
        if p.exists():
            return p

    # 3. Homebrew default path (Fallback)
    p = Path("/opt/homebrew/share/gromacs/top") / name
    if p.exists():
        return p

    raise FileNotFoundError(f"include file not found: {name}")


def read_atomtypes_recursive(top_file, visited=None):
    """Collect atomtypes by tracing all includes recursively"""
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

            # ---- Include processing ----
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
                        # Skip missing files like posre.itp with a warning (output to sys.stderr)
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

                # Remove inline comments (after ;)
                line_no_comment = line.split(';')[0].strip()
                if not line_no_comment:
                    continue

                parts = line_no_comment.split()
                if len(parts) >= 6:
                    atomtype = parts[0]
                    try:
                        # Force fields have different numbers of columns, so get the last two columns (sigma, epsilon)
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

                # Remove inline comments
                line_no_comment = line.split(';')[0].strip()
                if not line_no_comment:
                    continue

                if line_no_comment[0].isdigit():
                    parts = line_no_comment.split()
                    nr = int(parts[0])
                    type_name = parts[1]
                    # Get 'atom name' from the 5th column (index 4) of [ atoms ] as atom_label.
                    # Use type_name as a fallback if columns are missing.
                    atom_label = parts[4] if len(parts) > 4 else type_name
                    charge = float(parts[6])
                    atoms.append((nr, type_name, atom_label, charge))

    return atoms


def center_coords(coords):
    """Translate coordinates so the geometric center is at the origin"""
    n = len(coords)
    if n == 0:
        return coords
    cx = sum(c[0] for c in coords) / n
    cy = sum(c[1] for c in coords) / n
    cz = sum(c[2] for c in coords) / n
    return [(x - cx, y - cy, z - cz) for x, y, z in coords]


def write_pdb(filename, atoms, coords, n_atoms):
    """Write coordinates to a PDB file"""
    with open(filename, 'w') as f:
        for i in range(n_atoms):
            atom_label = atoms[i][2]
            x, y, z = coords[i]
            # Element symbol (guess from first letter, keep alphabetic only)
            element = ''.join([c for c in atom_label if c.isalpha()])[:1]
            # PDB ATOM format
            f.write(f"ATOM  {i+1:>5} {atom_label:<4} UNK A   1    {x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {element:>2}\n")
        f.write("END\n")


def main():
    parser = argparse.ArgumentParser(description="Convert GROMACS top/gro to RISMiCal prm.")
    parser.add_argument("-c", "--centering", action="store_true", help="Center the molecule geometric center to the origin and output a PDB file")
    parser.add_argument("top_file", help="Input topology file (.top)")
    parser.add_argument("gro_file", help="Input coordinate file (.gro)")
    
    args = parser.parse_args()

    top_file = args.top_file
    gro_file = args.gro_file

    atomtypes = read_atomtypes_recursive(top_file)
    atoms = read_atoms_recursive(top_file)
    n_atoms, coords = read_gro(gro_file)

    # Apply centering if requested
    if args.centering:
        coords = center_coords(coords)
        pdb_filename = f"{Path(gro_file).stem}_centering.pdb"
        write_pdb(pdb_filename, atoms, coords, n_atoms)
        # Output message to stderr so it doesn't corrupt redirected .prm output
        print(f"[Info] Centering applied. Centered PDB saved to: {pdb_filename}", file=sys.stderr)

    # Line 1: number_of_atoms name_of_molecule(optional)
    molecule_name = Path(top_file).stem
    print(f"{n_atoms} {molecule_name}")

    # Output for Line 2 and onwards
    for i, (nr, type_name, atom_label, charge) in enumerate(atoms[:n_atoms]):
        if type_name not in atomtypes:
            raise KeyError(f"atomtype '{type_name}' not found in any include")

        sigma, epsilon = atomtypes[type_name]
        x, y, z = coords[i]

        # atom_label LJ_sigma[Angs] LJ_epsilon[J/mol] point_charge[e] coordinate_x, y, z [Angs]
        print(f"{atom_label:<8s} {sigma:8.4f} {epsilon:10.4f} {charge:8.4f} {x:8.4f} {y:8.4f} {z:8.4f}")


if __name__ == "__main__":
    main()