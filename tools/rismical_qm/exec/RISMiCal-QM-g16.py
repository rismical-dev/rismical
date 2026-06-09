#!/usr/bin/env python3
"""
RISMiCal-QM-g16.py
A wrapper script for performing QM/3D-RISM calculations using RISMiCal and Gaussian 16.
[MODIFIED] Added Picard mixing (under-relaxation) for the electrostatic potential.
[CRITICAL FIX] Fixed Cube file writing format to properly respect z-axis line breaks.
[RESTORED] Restored all detailed comments and console output messages for usability.
[MODIFIED] Output units changed to kJ/mol. Energy equation restored to standard QM/MM scheme.
[NEW] Enhanced FC mode output to dynamically extract and display TD-DFT excited states.
[CRITICAL FIX] Robust target energy parsing and explicit coordinate writing for E(UU) extraction.
[CRITICAL FIX] Maintained NoSymm in E(UU) extraction to prevent Guess=Read orientation mismatch crashes.
[NEW] Added safe handling of .chk files during E(UU) extraction to prevent density contamination.
[NEW] Added log truncation (save_truncated_gout) to save disk space by discarding massive point-charge lists.
[MODIFIED] Removed false-positive termination warning for E(UU) extraction (MaxCycle=1 intentionally does not converge).
[NEW] Re-structured final summary output and added rigorous E_QM-MM extraction from polarized density.
"""

import os
import sys
import shutil
import subprocess
import re
import numpy as np
from scipy.spatial.distance import cdist

# ==============================================================================
# Commands and Paths Configuration
# ==============================================================================
G16_CMD      = "g16"
FORMCHK_CMD  = "formchk"
CUBEGEN_CMD  = "cubegen"
RISMICAL_CMD = "rismical.x"

# Physical constants for energy and length conversions (MODIFIED to kJ/mol)
HARTREE_TO_KJMOL = 2625.4996394799
COULOMB_TO_KJMOL = 1389.35456      # 1 e^2 / Angstrom -> kJ/mol
ANG_TO_BOHR      = 1.8897261246    # Angstrom -> Bohr

# Change Fortran double precision notation 'd' to 'e' for Python float conversion
def parse_fortran_float(val_str):
    try:
        return float(str(val_str).lower().replace('d', 'e'))
    except ValueError:
        return 0.0

# Parse the qmpart string to create a list of indices for QM atoms
def parse_qmpart(qmpart_str, total_atoms):
    # If qmpart is not given, all atoms are treated as QM
    if not qmpart_str or not str(qmpart_str).strip():
        return list(range(1, total_atoms + 1))
    
    indices = []
    parts = str(qmpart_str).split(',')
    last_idx = 0
    for p in parts:
        if not p.strip(): continue
        val = int(p.strip())
        if val > 0:
            indices.append(val)
            last_idx = val
        elif val < 0:
            # Handle negative values as a range from the last positive index
            indices.extend(range(last_idx + 1, abs(val) + 1))
            last_idx = abs(val)
    return sorted(list(set(indices)))

# Extract parameters from a specific namelist block string
def extract_namelist(block_str, params):
    placeholders = []
    def repl(m):
        placeholders.append(m.group(0))
        return f"__QUOTE_{len(placeholders)-1}__"
        
    masked_str = re.sub(r'"[^"]*"|\'[^\']*\'', repl, block_str)
    matches = list(re.finditer(r'([a-zA-Z_]\w*)\s*=', masked_str))
    
    for i, match in enumerate(matches):
        key = match.group(1).lower()
        start_idx = match.end()
        end_idx = matches[i+1].start() if i+1 < len(matches) else len(masked_str)
        val_str = masked_str[start_idx:end_idx].strip().rstrip(',')
        for j, p in enumerate(placeholders): 
            val_str = val_str.replace(f"__QUOTE_{j}__", p)
        if (val_str.startswith('"') and val_str.endswith('"')) or (val_str.startswith("'") and val_str.endswith("'")):
            val_str = val_str[1:-1]
        params[key] = val_str

# Read the original RISMiCal input file and separate namelists and atom data
def read_input_file(inp_file):
    params, udata, udata_line_indices = {}, [], []
    with open(inp_file, 'r', errors='replace') as f: lines = f.readlines()
    
    in_rismicalqm, in_grid3d, in_udata = False, False, False
    rismicalqm_block_str, grid3d_block_str = "", ""
    
    for i, line in enumerate(lines):
        clean_line = line.split('!')[0].strip()
        if not clean_line: continue
        upper_line = clean_line.upper()
        
        # Detect start of blocks
        if upper_line.startswith('$RISMICALQM') or upper_line.startswith('&RISMICALQM'):
            in_rismicalqm = True
            header_len = len('$RISMICALQM') if upper_line.startswith('$RISMICALQM') else len('&RISMICALQM')
            rismicalqm_block_str += clean_line[header_len:] + " "
            continue
        elif upper_line.startswith('$GRID3D') or upper_line.startswith('&GRID3D'):
            in_grid3d = True
            header_len = len('$GRID3D') if upper_line.startswith('$GRID3D') else len('&GRID3D')
            grid3d_block_str += clean_line[header_len:] + " "
            continue
        elif upper_line.startswith('$UDATA') or upper_line.startswith('&UDATA'):
            in_udata = True; continue
            
        # Detect end of blocks
        if in_rismicalqm and (upper_line == '$END' or upper_line == '/'):
            in_rismicalqm = False; continue
        elif in_grid3d and (upper_line == '$END' or upper_line == '/'):
            in_grid3d = False; continue
        elif in_udata and (upper_line == '$END' or upper_line == '/'):
            in_udata = False; continue
        
        # Append data to corresponding blocks
        if in_rismicalqm: rismicalqm_block_str += clean_line + " "
        elif in_grid3d: grid3d_block_str += clean_line + " "
        elif in_udata:
            parts = clean_line.split()
            if len(parts) >= 7: 
                udata.append(parts)
                udata_line_indices.append(i)
                
    extract_namelist(rismicalqm_block_str, params)
    extract_namelist(grid3d_block_str, params)
    return params, udata, udata_line_indices, lines

# Generate the Gaussian job file (.gjf) dynamically based on parameters
def write_gaussian_input(gjf_file, params, udata, qm_indices, iter_num, qv_file):
    qmopt_raw = params.get('qmopt', '')
    raw_lines = qmopt_raw.split('\\')
    link0_lines, route_line, title_line, chg_mult = [], "", "RISMiCal-QM Job", "0 1"
    
    for line in raw_lines:
        line = line.strip()
        if line.startswith('%'): link0_lines.append(line)
        elif line.startswith('#'): route_line = line
        elif line and line[0].isdigit() and ' ' in line: chg_mult = line
        elif line: title_line = line
            
    # Ensure MK population analysis is requested for ESP charges
    if not re.search(r'(?i)\bpop=', route_line): route_line += " Pop=MK"
            
    # iter_num == -1 specifies Pure QM Vacuum (No MM, No Charge)
    if iter_num == -1:
        route_line = re.sub(r'(?i)\bcharge\b', '', route_line)
        route_line = re.sub(r'(?i)\bnosymm\b', '', route_line)
    else:
        if not re.search(r'(?i)\bcharge\b', route_line): route_line += " Charge NoSymm"

    with open(gjf_file, 'w') as f:
        for l in link0_lines: f.write(f"{l}\n")
        f.write(f"{route_line}\n\n{title_line}\n\n{chg_mult}\n")
        
        # Write QM atoms coordinates
        for idx in qm_indices:
            atom = udata[idx - 1]
            x, y, z = parse_fortran_float(atom[4]), parse_fortran_float(atom[5]), parse_fortran_float(atom[6])
            f.write(f"{atom[0]:<2s}  {x:11.6f}  {y:11.6f}  {z:11.6f}\n")
        f.write("\n")
        
        # Link the external charge file if applicable
        if iter_num >= 0 and os.path.exists(qv_file): 
            f.write(f"@{qv_file}\n\n")

# [NEW] Helper function to parse all computed TD-DFT excited states
def extract_td_info(gout_file):
    """
    Parses the Gaussian output to extract the ground state SCF energy and 
    the excitation energy (eV) / absolute energy (Hartree) of all computed TD states.
    """
    scf_energy = 0.0
    td_states = []
    if not os.path.exists(gout_file): return scf_energy, td_states
    
    with open(gout_file, 'r', errors='replace') as f:
        lines = f.readlines()
        
    # Get ground state SCF energy first (used as baseline)
    for line in reversed(lines):
        if 'SCF Done:' in line:
            try:
                scf_energy = float(line.split()[4])
            except ValueError:
                pass
            break
            
    # Extract TD-DFT excited states dynamically
    for line in lines:
        m = re.search(r'Excited State\s+(\d+):\s+(.*?)\s+([\d\.]+)\s+eV', line)
        if m:
            st = int(m.group(1))
            sym = m.group(2).strip()
            ev = float(m.group(3))
            hartree = scf_energy + ev / 27.211386245988
            td_states.append({'state': st, 'symmetry': sym, 'ev': ev, 'hartree': hartree})
            
    return scf_energy, td_states

# [MODIFIED] Unified function to reliably extract total energy for Ground or Target Excited State
def read_target_energy(gout_file):
    """ Reads target state energy (SCF Done or TD Total Energy) in Hartree safely. """
    if not os.path.exists(gout_file): return 0.0
    with open(gout_file, 'r', errors='replace') as f:
        lines = f.readlines()
    
    # Check for TD energy formats
    for line in reversed(lines):
        if 'E(TD-HF/TD-DFT)' in line or 'Total Energy, E(' in line:
            try:
                val = line.split('=')[-1].split()[0]
                return float(val)
            except (IndexError, ValueError):
                pass
        
    # Fallback to ground state SCF energy
    for line in reversed(lines):
        if 'SCF Done:' in line:
            try:
                return float(line.split()[4])
            except (IndexError, ValueError):
                pass
    return 0.0

# [NEW] Safe Energy Extraction: General function for both E(UU) and E_QM-MM extraction
def extract_energy_from_density(base_name, params, chk_file, udata, qm_indices, qv_file="", suffix="euu"):
    """
    Evaluates the energy of the polarized density under specific background charges.
    Explicitly maintains NoSymm to prevent coordinate reorientation errors, 
    and uses an isolated .chk copy to prevent main density contamination.
    """
    job_gjf = f"{base_name}_{suffix}.gjf"
    job_gout = f"{base_name}_{suffix}.gout"
    job_chk = f"{base_name}_{suffix}.chk"
    
    # Safely copy the main checkpoint to avoid altering it during extraction
    if os.path.exists(chk_file):
        shutil.copy(chk_file, job_chk)
    
    qmopt_raw = params.get('qmopt', '')
    raw_lines = qmopt_raw.split('\\')
    link0_lines, route_line, title_line, chg_mult = [], "", f"Energy Extraction ({suffix})", "0 1"
    
    # Parse original inputs
    for line in raw_lines:
        line = line.strip()
        if line.startswith('%'): 
            if 'chk=' not in line.lower():
                link0_lines.append(line)
        elif line.startswith('#'): route_line = line
        elif line and line[0].isdigit() and ' ' in line: chg_mult = line
        elif line: title_line = line
            
    # Target the copied isolated checkpoint
    link0_lines.append(f"%chk={job_chk}")
    
    # Strip Charge keyword by default, maintaining NoSymm
    route_line = re.sub(r'(?i)\bcharge\b', '', route_line)
    
    # If qv_file is provided (e.g., MM only charges), reinstate Charge keyword
    if qv_file and os.path.exists(qv_file):
        if not re.search(r'(?i)\bcharge\b', route_line):
            route_line += " Charge"
    
    # Add Guess=Read and 1-step cycle limits
    if not re.search(r'(?i)\bguess=read\b', route_line): 
        route_line += " Guess=Read"
    route_line += " SCF=(MaxCycle=1)"
    
    # Construct GJF file explicitly with coordinates
    with open(job_gjf, 'w') as f:
        for l in link0_lines: f.write(f"{l}\n")
        f.write(f"{route_line}\n\n{title_line}\n\n{chg_mult}\n")
        
        # Write exact QM atoms coordinates securely
        for idx in qm_indices:
            atom = udata[idx - 1]
            x, y, z = parse_fortran_float(atom[4]), parse_fortran_float(atom[5]), parse_fortran_float(atom[6])
            f.write(f"{atom[0]:<2s}  {x:11.6f}  {y:11.6f}  {z:11.6f}\n")
        f.write("\n")
        
        if qv_file and os.path.exists(qv_file):
            f.write(f"@{qv_file}\n\n")
            
    subprocess.run([G16_CMD, job_gjf, job_gout])
    
    # Read the extracted energy
    e_hartree = read_target_energy(job_gout)
    return e_hartree * HARTREE_TO_KJMOL

# [NEW] Helper function to truncate large Gaussian logs
def save_truncated_gout(src_file, dest_file):
    """
    Saves a truncated version of the Gaussian output to save disk space.
    Only keeps the lines from the final 'SCF Done:' onwards, discarding 
    the massive external charge coordinate lists.
    """
    if not os.path.exists(src_file): return
    with open(src_file, 'r', errors='replace') as f:
        lines = f.readlines()
        
    start_idx = 0
    # Find the last occurrence of SCF Done to catch the final converged iteration
    for i in range(len(lines)-1, -1, -1):
        if 'SCF Done:' in lines[i]:
            start_idx = i
            break
            
    with open(dest_file, 'w', errors='replace') as f:
        f.writelines(lines[start_idx:])

# Write a .qv file containing ONLY MM atoms (used for QM+MM Gas Phase calculations)
def write_mm_qv(qv_file, udata, mm_indices):
    with open(qv_file, 'w') as f:
        for i in mm_indices:
            atom = udata[i-1]
            x, y, z = parse_fortran_float(atom[4]), parse_fortran_float(atom[5]), parse_fortran_float(atom[6])
            q = parse_fortran_float(atom[3])
            f.write(f" {x:11.6f} {y:11.6f} {z:11.6f}  {q:.6e}\n")

# Check if Gaussian terminated normally by reading the last few lines of the log
def check_gaussian_termination(gout_file):
    if not os.path.exists(gout_file): return False
    with open(gout_file, 'r', errors='replace') as f:
        for line in f.readlines()[-20:]:
            if "Normal termination" in line: return True
    return False

# Check if RISMiCal terminated normally by reading the last few lines of the output
def check_rismical_termination(rsmout_file):
    if not os.path.exists(rsmout_file): return False
    with open(rsmout_file, 'r', errors='replace') as f:
        for line in f.readlines()[-30:]:
            if "RISMiCal computation is completed normally" in line: return True
    return False

# Read RESP/MK charges from the Gaussian output file
def read_charges_from_gout(gout_file, natoms_qm):
    with open(gout_file, 'r', errors='replace') as f: lines = f.readlines()
    start_idx = -1
    
    for i in range(len(lines)-1, -1, -1):
        if "Charges from ESP fit" in lines[i] or "Fitting point charges" in lines[i]:
            start_idx = i; break
            
    # Fallback to Mulliken charges if ESP fit fails or is not available
    if start_idx == -1:
        for i in range(len(lines)-1, -1, -1):
            if "Mulliken charges:" in lines[i]: start_idx = i; break
            
    charges = []
    if start_idx != -1:
        for line in lines[start_idx:]:
            parts = line.split()
            if len(parts) == 3 and parts[0].isdigit() and not parts[1].isdigit():
                try: charges.append(float(parts[2]))
                except: pass
            if len(charges) == natoms_qm: break
            
    return charges if len(charges) == natoms_qm else [0.0]*natoms_qm

# Run Gaussian cubegen utility to generate the ASCII cube file for electrostatic potential
def run_cubegen(fchk_file, ascii_cube_file, ngrid3d, rdelta3d, params):
    # Define exact grid origin and dimensions to match RISMiCal internal grids
    origin_ang = -rdelta3d * (ngrid3d / 2.0)
    grid_spec = f"0 {origin_ang:.6f} {origin_ang:.6f} {origin_ang:.6f}\n"
    grid_spec += f"{ngrid3d} {rdelta3d:.6f} 0.0 0.0\n"
    grid_spec += f"{ngrid3d} 0.0 {rdelta3d:.6f} 0.0\n"
    grid_spec += f"{ngrid3d} 0.0 0.0 {rdelta3d:.6f}\n"
    
    potential_type = "potential=scf"
    # Support for excited state potential generation
    if re.search(r'(?i)\bdensity=current\b', params.get('qmopt', '')):
        potential_type = "potential=current"
    
    result = subprocess.run([CUBEGEN_CMD, "0", potential_type, fchk_file, ascii_cube_file, "-1"],
                            input=grid_spec, text=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
                            
    if not os.path.exists(ascii_cube_file):
        print(f"\n  [ERROR] {ascii_cube_file} was not generated by cubegen.")
        if result.stderr: print(f"  [cubegen STDERR] {result.stderr.strip()}")
        sys.exit(1)

# Add electrostatic potential from MM atoms to the generated QM cube, applying Picard mixing if needed
def add_mm_potential_to_cube(ascii_cube_file, udata, mm_indices, mix_alpha=1.0, prev_pot=None):
    """
    Reads the ASCII cube, adds MM potential, applies Picard mixing if requested,
    and overwrites the cube file accurately matching Fortran read rules.
    """
    with open(ascii_cube_file, 'r', errors='replace') as f:
        lines = f.readlines()

    natoms = int(lines[2].split()[0])
    abs_natoms = abs(natoms)
    header_end = 6 + abs_natoms

    # Extract origin and steps
    x0, y0, z0 = map(float, lines[2].split()[1:4])
    nx, dx = float(lines[3].split()[0]), float(lines[3].split()[1])
    ny, dy = float(lines[4].split()[0]), float(lines[4].split()[2])
    nz, dz = float(lines[5].split()[0]), float(lines[5].split()[3])
    nx, ny, nz = int(nx), int(ny), int(nz)

    # Read original QM potential (flat array)
    data_flat = []
    for line in lines[header_end:]:
        data_flat.extend(map(float, line.split()))
    qm_pot = np.array(data_flat, dtype=np.float64)

    # Generate Grid Vectors and calculate MM Potential V = sum(q/r)
    v_mm = np.zeros_like(qm_pot)
    if mm_indices:
        mm_c, mm_q = [], []
        for i in mm_indices:
            atom = udata[i-1]
            mm_c.append([parse_fortran_float(atom[4]) * ANG_TO_BOHR, 
                         parse_fortran_float(atom[5]) * ANG_TO_BOHR, 
                         parse_fortran_float(atom[6]) * ANG_TO_BOHR])
            mm_q.append(parse_fortran_float(atom[3]))

        X = x0 + np.arange(nx) * dx
        Y = y0 + np.arange(ny) * dy
        Z = z0 + np.arange(nz) * dz
        XX, YY, ZZ = np.meshgrid(X, Y, Z, indexing='ij')
        XX_flat, YY_flat, ZZ_flat = XX.ravel(), YY.ravel(), ZZ.ravel()

        for c, q in zip(mm_c, mm_q):
            d_sq = (XX_flat - c[0])**2 + (YY_flat - c[1])**2 + (ZZ_flat - c[2])**2
            d_sq = np.maximum(d_sq, 1e-12) # Prevent div by zero if MM atom lies on grid
            v_mm += q / np.sqrt(d_sq)

    # New complete potential (QM + MM)
    new_pot = qm_pot + v_mm

    # Apply Picard Mixing: P_mix = alpha * P_new + (1 - alpha) * P_prev
    mixed_pot = new_pot
    if prev_pot is not None and mix_alpha < 1.0:
        print(f"  [INFO] Applying Picard mixing with alpha = {mix_alpha}")
        mixed_pot = mix_alpha * new_pot + (1.0 - mix_alpha) * prev_pot

    # Helper function to segment array for Fortran 6E13.5 formatting
    def chunker(seq, size):
        return (seq[pos:pos + size] for pos in range(0, len(seq), size))

    # Reshape to 3D to ensure line breaks occur at the end of each z-axis loop
    mixed_pot_3d = mixed_pot.reshape((nx, ny, nz))

    with open(ascii_cube_file, 'w') as f:
        f.writelines(lines[:header_end])
        for ix in range(nx):
            for iy in range(ny):
                row = mixed_pot_3d[ix, iy, :]
                for chunk in chunker(row, 6):
                    line = "".join([f"{v:13.5E}" for v in chunk]) + "\n"
                    f.write(line)
            
    return mixed_pot

# Update RISMiCal input file with new charges for QM atoms
def update_rismical_input(inp_file, lines, udata, udata_line_indices, qm_indices, charges):
    for c_idx, idx in enumerate(qm_indices):
        row = udata[idx-1]
        p1, p2 = parse_fortran_float(row[1]), parse_fortran_float(row[2])
        x, y, z = parse_fortran_float(row[4]), parse_fortran_float(row[5]), parse_fortran_float(row[6])
        lines[udata_line_indices[idx-1]] = f"{row[0]:<7s}{p1:8.4f}{p2:10.4f}{charges[c_idx]:12.6f}     {x:8.4f}   {y:8.4f}   {z:8.4f}\n"
    with open(inp_file, 'w', errors='replace') as f: f.writelines(lines)

# Filter solvent charges to avoid Coulomb singularity and calculate MM-Solvent interaction energy
def process_qv_and_calc_emv(qv_file, udata, qvcutoff, qvcore, mm_indices):
    qv_coords, qv_charges = [], []
    
    # Read solvent charges from .qv
    with open(qv_file, 'r', errors='replace') as fin:
        for line in fin:
            clean = line.split('!')[0].strip() 
            if not clean: continue
            parts = clean.split()
            if len(parts) >= 4:
                try:
                    q = float(parts[3])
                    if abs(q) >= qvcutoff:
                        qv_coords.append([float(parts[0]), float(parts[1]), float(parts[2])])
                        qv_charges.append(q)
                except ValueError: continue
                    
    if not qv_coords: return 0.0
        
    qv_c_arr, qv_q_arr = np.array(qv_coords), np.array(qv_charges)
    all_solute_coords = [[parse_fortran_float(row[4]), parse_fortran_float(row[5]), parse_fortran_float(row[6])] for row in udata]
    
    # Apply qvcore distance filter
    dists = cdist(qv_c_arr, np.array(all_solute_coords))
    valid_mask = np.min(dists, axis=1) >= qvcore
    qv_c_valid, qv_q_valid = qv_c_arr[valid_mask], qv_q_arr[valid_mask]
    
    # Calculate electrostatic interaction energy between MM atoms and solvent
    e_mv_kjmol = 0.0
    if mm_indices and len(qv_c_valid) > 0:
        mm_c = np.array([[parse_fortran_float(udata[i-1][4]), parse_fortran_float(udata[i-1][5]), parse_fortran_float(udata[i-1][6])] for i in mm_indices])
        mm_q = np.array([parse_fortran_float(udata[i-1][3]) for i in mm_indices])
        inv_dists = np.where(cdist(mm_c, qv_c_valid) > 1e-6, 1.0/cdist(mm_c, qv_c_valid), 0.0)
        e_mv_kjmol = np.dot(mm_q, np.dot(inv_dists, qv_q_valid)) * COULOMB_TO_KJMOL
        
    # Write back the filtered solvent charges and MM charges to .qv file for Gaussian
    with open(qv_file, 'w') as fout:
        for c, q in zip(qv_c_valid, qv_q_valid):
            fout.write(f" {c[0]:11.6f} {c[1]:11.6f} {c[2]:11.6f}  {q:.6e}\n")
        if mm_indices:
            for i in mm_indices:
                atom = udata[i-1]
                x, y, z = parse_fortran_float(atom[4]), parse_fortran_float(atom[5]), parse_fortran_float(atom[6])
                q = parse_fortran_float(atom[3])
                fout.write(f" {x:11.6f} {y:11.6f} {z:11.6f}  {q:.6e}\n")
                
    return e_mv_kjmol

# Read RISMiCal summary output to get solvent free energy components (converted to kJ/mol)
def read_xmu(xmu_file):
    sfe_sc, se_es = 0.0, 0.0
    with open(xmu_file, 'r', errors='replace') as f:    
        for line in f:
            if 'SFE_SC=' in line: sfe_sc = float(line.split('=')[1].split('!')[0].strip()) / 1000.0
            elif 'SE_ES=' in line: se_es = float(line.split('=')[1].split('!')[0].strip()) / 1000.0
    return sfe_sc, se_es

# Extract total energy from Gaussian output log using robust parsing
def read_eqm(gout_file):
    return read_target_energy(gout_file) * HARTREE_TO_KJMOL

# ==============================================================================
# Main Execution Block
# ==============================================================================
def main():
    # Detect FC mode
    fc_mode = False
    if "-FC" in sys.argv:
        fc_mode = True
        sys.argv.remove("-FC")

    if len(sys.argv) < 2: 
        print("Usage: python RISMiCal-QM-g16.py <input_file> [-FC]"); sys.exit(1)
        
    inp_file = sys.argv[1]
    base_name = os.path.splitext(inp_file)[0]
    
    # Define file paths
    gjf_file, chk_file, fchk_file = f"{base_name}.gjf", f"{base_name}.chk", f"{base_name}.fchk"
    gout_file, ascii_cube = f"{base_name}.gout", f"{base_name}.cub"
    xmu_file, qv_file, org_inp_file = f"{base_name}.xmu", f"{base_name}.qv", f"{base_name}.org.inp"
    
    # Backup original input
    if not os.path.exists(org_inp_file): shutil.copy(inp_file, org_inp_file)
    
    params, udata, udata_indices, all_lines = read_input_file(inp_file)
    
    # Ensure a global checkpoint file is declared for all jobs in the sequence
    if '%chk' not in params.get('qmopt', '').lower():
        params['qmopt'] = f"%chk={chk_file}\\" + params.get('qmopt', '')
        
    # Read calculation configuration
    scfconv   = parse_fortran_float(params.get('scfconv', '1e-4'))
    qvcutoff  = parse_fortran_float(params.get('qvcutoff', '1e-6'))
    # Set qvcore in the input file (e.g., qvcore = 0.1) if testing standard method
    qvcore    = parse_fortran_float(params.get('qvcore', '0.5'))
    ngrid3d   = int(parse_fortran_float(params.get('ngrid3d', '128')))
    rdelta3d  = parse_fortran_float(params.get('rdelta3d', '0.5'))
    
    # Extract mix_alpha from input params for Picard method (default is 1.0, no mixing)
    mix_alpha = parse_fortran_float(params.get('mix_alpha', '1.0'))
    
    total_atoms = len(udata)
    qm_indices = parse_qmpart(params.get('qmpart', ''), total_atoms)
    mm_indices = [i for i in range(1, total_atoms + 1) if i not in qm_indices]
    
    print(f"--- RISMiCal-QM-g16 Started ({base_name}) ---")
    print(f" Total Atoms: {total_atoms} (QM: {len(qm_indices)}, MM: {len(mm_indices)})")
    print(f" Distance Filter (qvcore): {qvcore} Angstrom")
    if mix_alpha < 1.0:
        print(f" Convergence Acceleration: Picard Mixing (alpha = {mix_alpha})")

    # =========================================================
    # Franck-Condon (FC) Mode Execution
    # =========================================================
    if fc_mode:
        print("\n--- Franck-Condon (FC) State Calculation ---")
        if not os.path.exists(qv_file):
            print(f"  [ERROR] FC mode requires an existing '{qv_file}'.")
            sys.exit(1)
        print(f"  Using existing frozen solvent/MM charges from: {qv_file}")
            
        write_gaussian_input(gjf_file, params, udata, qm_indices, iter_num=1, qv_file=qv_file)
        
        print("  Running Gaussian 16...")
        subprocess.run([G16_CMD, gjf_file, gout_file])
        if not check_gaussian_termination(gout_file):
            print(f"  [ERROR] Gaussian 16 terminated abnormally. Check {gout_file}"); sys.exit(1)
            
        print("  >>> FC Calculation Completed Successfully! <<<")
        
        # Output energy and excitation energy (eV) for each state
        scf_e, td_states = extract_td_info(gout_file)
        if td_states:
            print("\n  [Excited States Information]")
            print(f"  Ground State (SCF) Energy : {scf_e:.6f} Hartree")
            for st in td_states:
                print(f"  State {st['state']:2d} ({st['symmetry']:<12s}) : {st['hartree']:.6f} Hartree  ({st['ev']:.4f} eV)")
                
        # Extract pure internal E(UU)
        e_uu_kjmol = extract_energy_from_density(base_name, params, chk_file, udata, qm_indices, qv_file="", suffix="euu")
        
        # Extract QM-MM interaction energy under solvent polarized density
        e_qm_mm_kjmol = 0.0
        if mm_indices:
            mm_qv_file = f"{base_name}_mm.qv"
            write_mm_qv(mm_qv_file, udata, mm_indices)
            e_qmmm_solv_kjmol = extract_energy_from_density(base_name, params, chk_file, udata, qm_indices, qv_file=mm_qv_file, suffix="eqmmm")
            e_qm_mm_kjmol = e_qmmm_solv_kjmol - e_uu_kjmol
            
        e1_hartree = read_target_energy(gout_file)
        e1_kjmol = e1_hartree * HARTREE_TO_KJMOL
        
        # QM-Solvent interaction is the residual after removing internal and QM-MM energies
        e_qm_solv_kjmol = e1_kjmol - e_uu_kjmol - e_qm_mm_kjmol
        
        print("\n  [Target State (Root) Properties]")
        print(f"  QM Internal Energy (E_UU)    : {e_uu_kjmol:.4f} kJ/mol")
        if mm_indices:
            print(f"  QM-MM Int. (E_QM-MM)         : {e_qm_mm_kjmol:.4f} kJ/mol")
        print(f"  QM-Solvent Int. (E_QM-v)     : {e_qm_solv_kjmol:.4f} kJ/mol")
        
        # Output effective charges of each atom for the target state
        fc_chg = read_charges_from_gout(gout_file, len(qm_indices))
        print("\n  Effective Charges (ESP/MK) for Target State:")
        for i, q in zip(qm_indices, fc_chg):
            atom_name = udata[i-1][0]
            print(f"    Atom {i:2d} ({atom_name:2s}) : {q:9.6f} e")

        print(f"\n  Check '{base_name}.gout.fc' for excited state properties.")
        save_truncated_gout(gout_file, f"{base_name}.gout.fc")
        sys.exit(0)

    # ---------------------------------------------------------
    # Pre-Step 1: Pure QM Vacuum
    # ---------------------------------------------------------
    print("\n--- Pre-Step 1: QM Vacuum (E_gas) ---")
    write_gaussian_input(gjf_file, params, udata, qm_indices, iter_num=-1, qv_file="")
    subprocess.run([G16_CMD, gjf_file, gout_file])
    if not check_gaussian_termination(gout_file):
        print("  [ERROR] Gaussian 16 terminated abnormally in Pre-Step 1."); sys.exit(1)
    e_gas = read_eqm(gout_file)
    save_truncated_gout(gout_file, f"{base_name}.gout.gas")

    # ---------------------------------------------------------
    # Pre-Step 2: QM + MM Vacuum
    # ---------------------------------------------------------
    print("--- Pre-Step 2: QM+MM Vacuum (E_QMMM_gas) ---")
    e_qmmm_gas = e_gas
    e_qmmm_int = 0.0
    if mm_indices:
        write_mm_qv(qv_file, udata, mm_indices)
        write_gaussian_input(gjf_file, params, udata, qm_indices, iter_num=0, qv_file=qv_file)
        subprocess.run([G16_CMD, gjf_file, gout_file])
        if not check_gaussian_termination(gout_file):
            print("  [ERROR] Gaussian 16 terminated abnormally in Pre-Step 2."); sys.exit(1)
        e_qmmm_gas = read_eqm(gout_file)
        e_qmmm_int = e_qmmm_gas - e_gas
        save_truncated_gout(gout_file, f"{base_name}.gout.qmmm_gas")
    else:
        print("  No MM atoms specified. Skipping...")

    # ---------------------------------------------------------
    # Initializing Solvent Distribution (Iteration 0 for RISM)
    # ---------------------------------------------------------
    print("\n--- Initializing Solvent Distribution ---")
    new_chg = read_charges_from_gout(gout_file, len(qm_indices))
    subprocess.run([FORMCHK_CMD, chk_file, fchk_file], stdout=subprocess.DEVNULL)
    run_cubegen(fchk_file, ascii_cube, ngrid3d, rdelta3d, params)
    
    # Save returned potential to prev_pot to use in next iteration
    prev_pot = add_mm_potential_to_cube(ascii_cube, udata, mm_indices, mix_alpha=1.0, prev_pot=None)
    
    update_rismical_input(inp_file, all_lines, udata, udata_indices, qm_indices, new_chg)
    current_rsmout = f"{base_name}.rsmout.0"
    with open(current_rsmout, "w", errors='replace') as f_rsm: 
        subprocess.run([RISMICAL_CMD, "3d", inp_file], stdout=f_rsm, stderr=subprocess.STDOUT)
    if not check_rismical_termination(current_rsmout):
        print(f"  [ERROR] RISMiCal did not terminate normally. Check {current_rsmout}"); sys.exit(1)

    # ---------------------------------------------------------
    # Main SCF Loop
    # ---------------------------------------------------------
    prev_e = 0.0
    iter_num = 1
    print("\n--- Starting QM/MM/3D-RISM SCF Loop ---")
    while True:
        print(f"\n[Iteration {iter_num}]")
        
        sfe, se = read_xmu(xmu_file)
        
        # e_mv (MM-Solvent interaction) is calculated and charges are filtered, 
        # but e_mv itself is no longer displayed as per the user's specification.
        _ = process_qv_and_calc_emv(qv_file, udata, qvcutoff, qvcore, mm_indices)
        
        write_gaussian_input(gjf_file, params, udata, qm_indices, iter_num, qv_file)
        subprocess.run([G16_CMD, gjf_file, gout_file])
        if not check_gaussian_termination(gout_file):
            print(f"  [ERROR] Gaussian 16 terminated abnormally. Check {gout_file}"); sys.exit(1)
            
        # Extract pure internal E(UU) safely
        e_uu_kjmol = extract_energy_from_density(base_name, params, chk_file, udata, qm_indices, qv_file="", suffix="euu")
        
        # Extract QM-MM interaction energy under solvent polarized density
        e_qm_mm_kjmol = 0.0
        if mm_indices:
            mm_qv_file = f"{base_name}_mm.qv"
            write_mm_qv(mm_qv_file, udata, mm_indices)
            e_qmmm_solv_kjmol = extract_energy_from_density(base_name, params, chk_file, udata, qm_indices, qv_file=mm_qv_file, suffix="eqmmm")
            e_qm_mm_kjmol = e_qmmm_solv_kjmol - e_uu_kjmol
            
        # Save truncated log to avoid massive file sizes from point-charges
        save_truncated_gout(gout_file, f"{base_name}.gout.{iter_num}")
        
        new_chg = read_charges_from_gout(gout_file, len(qm_indices))
        
        # Calculate Total Free Energy: E_UU + E_QM-MM + SFE
        e_tot = e_uu_kjmol + e_qm_mm_kjmol + sfe
        
        print(f"  E(QM) Internal  = {e_uu_kjmol:.4f} kJ/mol")
        if mm_indices: print(f"  E_QM-MM         = {e_qm_mm_kjmol:.4f} kJ/mol")
        print(f"  E_QMMM-v        = {se:.4f} kJ/mol")
        print(f"  SFE_SC          = {sfe:.4f} kJ/mol")
        print(f"  G_TOTAL         = {e_tot:.4f} kJ/mol")
        
        d_e_h = abs(e_tot - prev_e) / HARTREE_TO_KJMOL
        print(f"  Delta E = {d_e_h:.6e} Hartree (Threshold: {scfconv})")
        
        if iter_num > 1 and d_e_h <= scfconv: break
        prev_e = e_tot
        
        # Potentials generation and mixing
        subprocess.run([FORMCHK_CMD, chk_file, fchk_file], stdout=subprocess.DEVNULL)
        run_cubegen(fchk_file, ascii_cube, ngrid3d, rdelta3d, params)
        
        # Pass mix_alpha and prev_pot for Picard mixing
        prev_pot = add_mm_potential_to_cube(ascii_cube, udata, mm_indices, mix_alpha=mix_alpha, prev_pot=prev_pot)
        
        update_rismical_input(inp_file, all_lines, udata, udata_indices, qm_indices, new_chg)
        current_rsmout = f"{base_name}.rsmout.{iter_num}"
        with open(current_rsmout, "w", errors='replace') as f_rsm: 
            subprocess.run([RISMICAL_CMD, "3d", inp_file], stdout=f_rsm, stderr=subprocess.STDOUT)
            
        if not check_rismical_termination(current_rsmout):
            print(f"  [ERROR] RISMiCal did not terminate normally in Iteration {iter_num}. Check {current_rsmout}"); sys.exit(1)
            
        iter_num += 1
        
    print("\n" + "="*55)
    print(" SUMMARY OF RISMiCal-QM CALCULATION")
    print("="*55)
    print(" >>> SCF Converged! <<<")
    print(f" QM Gas Energy (E_QM_gas)    : {e_gas:15.5f} kJ/mol")
    if mm_indices:
        print(f" QMMM Gas Energy (E_QMMM_gas): {e_qmmm_gas:15.5f} kJ/mol")
        print(f" QM-MM Interaction (gas)     : {e_qmmm_int:15.5f} kJ/mol")
    print("-" * 55)
    print(f" Total Free Energy (G_TOTAL)     : {e_tot:15.5f} kJ/mol")
    print(f" QM Internal Energy (E_UU)       : {e_uu_kjmol:15.5f} kJ/mol")
    if mm_indices:
        print(f" QM-MM Int. (E_QM-MM)            : {e_qm_mm_kjmol:15.5f} kJ/mol")
    print(f" QM/MM-Solvent ES Int. (E_QMMM-v): {se:15.5f} kJ/mol")
    print(f" Solvation Free E. (SFE_SC)      : {sfe:15.5f} kJ/mol")
    print("-" * 55)
    print(" Effective Charges (ESP/MK) for QM region:")
    for i, q in zip(qm_indices, new_chg):
        atom_name = udata[i-1][0]
        print(f"   Atom {i:2d} ({atom_name:2s}) : {q:9.6f} e")
    print("-" * 55)
    print(" Generated Files:")
    print(f"   Input Backup  : {org_inp_file}")
    print(f"   Gaussian INP  : {gjf_file}")
    print(f"   Gaussian LOGs : {base_name}.gout.* (Truncated)")
    print(f"   RISMiCal LOGs : {base_name}.rsmout.*")
    print(f"   ASCII Cube    : {ascii_cube}")
    print(f"   RISMiCal OUT  : {xmu_file}")
    print(f"   Ext. Charges  : {qv_file}")
    print("="*55)

if __name__ == "__main__": 
    main()
