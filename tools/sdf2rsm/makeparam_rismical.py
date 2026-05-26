import sys
import argparse
import traceback
from pathlib import Path

from openff.toolkit import Molecule, ForceField
from openff.units import unit

def main():
    parser = argparse.ArgumentParser(description="Extract MD parameters from an SDF file and output as RISMiCal text format.")
    parser.add_argument("input_file", help="Path to the input SDF file")
    parser.add_argument("mode", nargs='?', default="aa", help="'ua' for pseudo-United Atom mode (default is 'aa': All Atom)")
    args = parser.parse_args()

    input_path = Path(args.input_file)
    is_ua = (args.mode.lower() == "ua")

    if input_path.suffix.lower() != ".sdf":
        print(f"Error: Input file must be an SDF format (.sdf).", file=sys.stderr)
        sys.exit(1)
        
    if not input_path.exists():
        print(f"Error: File not found: {input_path.name}", file=sys.stderr)
        sys.exit(1)

    try:
        mode_name = 'United Atom' if is_ua else 'All Atom'
        print(f"Processing {input_path.name}... (Mode: {mode_name}, calculating AM1-BCC charges)")
        
        # 1. Load molecule and calculate charges
        molecule = Molecule.from_file(str(input_path))
        molecule.assign_partial_charges('am1bcc')
        charges = [q.m_as(unit.elementary_charge) for q in molecule.partial_charges]
        
        # 2. Load OpenFF Sage force field
        forcefield = ForceField("openff-2.1.0.offxml")
        labels = forcefield.label_molecules(molecule.to_topology())[0]
        vdw_labels = labels['vdW']
        
        # 3. Get coordinates
        if molecule.conformers is None or len(molecule.conformers) == 0:
            positions = [[0.0, 0.0, 0.0] for _ in range(molecule.n_atoms)]
        else:
            positions = molecule.conformers[0].m_as(unit.angstrom)

        # 4. Handle United Atom (UA) mode (merge aliphatic hydrogens into carbons)
        remove_indices = set()
        if is_ua:
            for atom in molecule.atoms:
                if atom.atomic_number == 6:  # Carbon
                    c_idx = atom.molecule_atom_index
                    for neighbor in atom.bonded_atoms:
                        if neighbor.atomic_number == 1:  # Hydrogen
                            h_idx = neighbor.molecule_atom_index
                            # Add H charge to C
                            charges[c_idx] += charges[h_idx]
                            # Mark H for removal
                            remove_indices.add(h_idx)
            print(f"Info: Merged {len(remove_indices)} hydrogens into their parent carbons.")

        # 5. Charge rounding error correction (ensure total charge is exactly an integer)
        # List of active atom indices (excluding removed hydrogens)
        active_indices = [i for i in range(molecule.n_atoms) if i not in remove_indices]
        
        # Target total charge (integer: 0, 1, -1, etc.)
        target_total_charge = round(sum(charges[i] for i in active_indices))
        
        # Round all charges to 4 decimal places initially
        for i in active_indices:
            charges[i] = round(charges[i], 4)
            
        # Calculate the difference caused by rounding
        current_total = sum(charges[i] for i in active_indices)
        diff = round(target_total_charge - current_total, 4)
        
        # If there is a difference, absorb it into the atom with the largest absolute charge
        if diff != 0.0:
            max_q_idx = max(active_indices, key=lambda i: abs(charges[i]))
            charges[max_q_idx] = round(charges[max_q_idx] + diff, 4)
            
            # Notify which atom (label) was adjusted
            max_q_label = molecule.atoms[max_q_idx].name or f"{molecule.atoms[max_q_idx].symbol}{max_q_idx+1}"
            print(f"Info: Adjusted total charge difference ({diff:+.4f} e) on atom '{max_q_label}'.")

        # 6. Collect output data
        atoms_data = []
        for i, atom in enumerate(molecule.atoms):
            # Skip hydrogens if UA mode is active
            if i in remove_indices:
                continue

            label = atom.name if atom.name else f"{atom.symbol}{i+1}"
            x, y, z = positions[i]
            
            # Get Lennard-Jones parameters
            vdw_param = vdw_labels[(i,)]
            sigma_angstrom = vdw_param.sigma.m_as(unit.angstrom)
            epsilon_j_mol = vdw_param.epsilon.m_as(unit.kilojoule_per_mole) * 1000.0
            q = charges[i]

            atoms_data.append((label, sigma_angstrom, epsilon_j_mol, q, x, y, z))

        # 7. Build RISMiCal format text
        molecule_name = input_path.stem + ("_ua" if is_ua else "")
        output_lines = [f"{len(atoms_data)}    {molecule_name}"]
        
        for data in atoms_data:
            label, sigma, epsilon, q, x, y, z = data
            
            # Format with Fortran 'd0' notation. 
            # q is already strictly rounded/corrected, so output as .4f directly.
            sigma_str = f"{sigma:.4f}d0"
            epsilon_str = f"{epsilon:.3f}d0"
            q_str = f"{q:.4f}d0"
            
            # Align columns for readability
            line = f"{label:<5} {sigma_str:>12} {epsilon_str:>12} {q_str:>12} {x:>14.7f} {y:>14.7f} {z:>14.7f}"
            output_lines.append(line)

        # 8. Save to file
        output_path = input_path.with_suffix(".txt")
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("\n".join(output_lines) + "\n")

        print(f"Success: RISMiCal parameters saved to {output_path.name}")

    except Exception as e:
        print(f"\nAn error occurred during processing:", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()