#!/usr/bin/env python3
"""
rismical-bind.py

This script automates the calculation of ligand binding free energies using the 3D-RISM theory,
mimicking an MM-PBSA workflow. It integrates ParmEd for topology/parameter parsing, 
MDAnalysis for trajectory handling and alignment, and RISMiCal for 3D-RISM calculations.

Workflow summary:
1. Parse RISMiCal input file to extract the $rismicalbind namelist.
2. Load topology parameters (LJ sigma/epsilon, charges) via ParmEd.
3. Load trajectory via MDAnalysis.
4. Iterate over trajectory frames:
   a. Perform RMSD fitting and/or centering if requested.
   b. If newtraj=.true., write the aligned coordinates to a new trajectory file.
   c. Calculate direct MM interaction energy (E_MM) between host and guest.
   d. Generate $UDATA blocks and input files for Host, Guest, and Complex.
      *If centereach=.true., dynamically shift the center of mass of EACH specific state 
       to (0,0,0) only when writing the .inp file to prevent box-edge artifacts.
   e. Execute RISMiCal externally for each state.
   f. Extract Solvation Free Energy (SFE_SC) and Correction_Term from .xmu outputs.
5. Compute the binding free energy (dG_bind) and output statistical summaries.
"""

import parmed as pmd
import MDAnalysis as mda
from MDAnalysis.analysis import align
from scipy.spatial.distance import cdist
import numpy as np
import subprocess
import re
import sys
import os
import csv

# =============================================================================
# Helper Functions for Input Parsing
# =============================================================================

def parse_atom_selection(sel_str):
    """
    Parses a custom atom selection string into a sorted list of 0-indexed integers.
    """
    if str(sel_str).strip() == '0' or not str(sel_str).strip():
        return []
    
    parts = str(sel_str).replace(" ", "").split(',')
    indices = []
    
    for i, p in enumerate(parts):
        val = int(p)
        if val < 0:
            start = int(parts[i-1])
            indices.extend(range(start + 1, abs(val) + 1))
        else:
            indices.append(val)
            
    return sorted(list(set(indices)))

def parse_rismicalbind(inp_text):
    """
    Extracts and parses the $rismicalbind namelist from the RISMiCal input file using regex.
    """
    params = {
        'host': [],
        'guest': [],
        'traj': [1, 'last', 1],
        'rmsfit': 'none',
        'centering': False,
        'newtraj': False,
        'centereach': False  # New option for independent state centering
    }
    
    match = re.search(r'\$rismicalbind(.*?)\$end', inp_text, re.IGNORECASE | re.DOTALL)
    if not match:
        raise ValueError("The $rismicalbind namelist was not found in the provided input file.")
        
    block = match.group(1)
    
    for line in block.split('\n'):
        line = line.strip()
        if not line or line.startswith('!'): 
            continue
        if '=' not in line: 
            continue
        
        key, val = line.split('=', 1)
        key = key.strip().lower()
        val = val.strip().replace('"', '').replace("'", "")
        
        if key == 'host':
            params['host'] = parse_atom_selection(val)
        elif key == 'guest':
            params['guest'] = parse_atom_selection(val)
        elif key == 'traj':
            parts = val.replace(' ', '').split(',')
            params['traj'][0] = int(parts[0])
            params['traj'][1] = parts[1] if parts[1].lower() == 'last' else int(parts[1])
            if len(parts) > 2:
                params['traj'][2] = int(parts[2])
        elif key == 'rmsfit':
            params['rmsfit'] = val.lower()
        elif key == 'centering':
            params['centering'] = val.lower() in ['.true.', 'true', 't', '1']
        elif key == 'newtraj':
            params['newtraj'] = val.lower() in ['.true.', 'true', 't', '1']
        elif key == 'centereach':
            params['centereach'] = val.lower() in ['.true.', 'true', 't', '1']
            
    return params

def get_mda_selection(idx_list):
    if not idx_list: 
        return "none"
    return "bynum " + " ".join(map(str, idx_list))

# =============================================================================
# Energy Calculation & File I/O Functions
# =============================================================================

def calc_emm(host_pos, guest_pos, host_sig, guest_sig, host_eps, guest_eps, host_q, guest_q):
    if len(host_pos) == 0 or len(guest_pos) == 0:
        return 0.0
    
    dists = cdist(host_pos, guest_pos)
    dists[dists == 0.0] = np.inf 
    
    sig_ij = (host_sig[:, None] + guest_sig[None, :]) / 2.0
    eps_ij = np.sqrt(host_eps[:, None] * guest_eps[None, :])
    
    sr = sig_ij / dists
    e_lj = np.sum(4.0 * eps_ij * (sr**12 - sr**6))
    
    COULOMB_CONST = 1389354.57 
    q_ij = host_q[:, None] * guest_q[None, :]
    e_coulomb = np.sum(COULOMB_CONST * q_ij / dists)
    
    return e_lj + e_coulomb

def write_rism_inp(filename, base_inp_text, u, indices_0, names, sigmas, epsilons, charges, label, centereach=False):
    """
    Generates a RISMiCal input file by appending the $UDATA block to the base input text.
    If centereach is True, it independently shifts the COM of the requested atoms to (0,0,0)
    to prevent box-edge errors during the 3D-RISM calculation.
    """
    with open(filename, 'w') as f:
        f.write(base_inp_text)
        if not base_inp_text.endswith('\n'):
            f.write('\n')
            
        f.write(' $UDATA\n')
        f.write(f'{len(indices_0)}  {label}\n')
        
        # Extract coordinates only for the specific state atoms
        state_atoms = u.atoms[indices_0]
        pos = state_atoms.positions.copy()  # Copy to avoid modifying the MDAnalysis universe
        
        # Shift COM to origin if requested
        if centereach and len(indices_0) > 0:
            com = state_atoms.center_of_mass()
            pos = pos - com
            
        # Write coordinates
        for idx, i in enumerate(indices_0):
            f.write(f'{names[i]:<7} {sigmas[i]:8.4f} {epsilons[i]:10.4f} {charges[i]:9.4f} '
                    f'{pos[idx][0]:10.7f} {pos[idx][1]:10.7f} {pos[idx][2]:10.7f}\n')

        f.write(' $END\n')

def extract_xmu(filepath):
    sfe_sc, corr = 0.0, 0.0
    if not os.path.exists(filepath):
        print(f"Warning: Output file {filepath} not found. Returning zeros.")
        return sfe_sc, corr
        
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('SFE_SC='):
                sfe_sc = float(line.split('=')[1].split('!')[0].strip())
            elif line.startswith('Correction_Term='):
                corr = float(line.split('=')[1].split('!')[0].strip())
                
    return sfe_sc, corr

# =============================================================================
# Main Workflow
# =============================================================================

def main(top_file, traj_file, inp_file):
    # --- Step 1: Parse the RISMiCal input file ---
    print("Reading and parsing RISMiCal base input file...")
    with open(inp_file, 'r') as f:
        base_inp = f.read()
    
    params = parse_rismicalbind(base_inp)
    host_idx_0 = np.array(params['host']) - 1
    guest_idx_0 = np.array(params['guest']) - 1
    
    if len(guest_idx_0) > 0:
        complex_idx_0 = np.sort(np.concatenate((host_idx_0, guest_idx_0)))
    else:
        complex_idx_0 = host_idx_0
    
    # --- Step 2: Load and Extract Parameters via ParmEd ---
    print(f"Loading topology parameters from {top_file} via ParmEd...")
    top = pmd.load_file(top_file)
    
    names = np.array([a.name for a in top.atoms])
    sigmas = np.array([a.sigma for a in top.atoms])
    epsilons = np.array([a.epsilon for a in top.atoms]) * 4184.0 
    charges = np.array([a.charge for a in top.atoms])
    
    # --- Step 3: Setup MDAnalysis Trajectory and Slicing ---
    print(f"Loading trajectory {traj_file} via MDAnalysis...")
    
    forced_ncdf = False
    try:
        u = mda.Universe(top_file, traj_file)
        ref_u = mda.Universe(top_file, traj_file)
    except ValueError:
        print("  -> GROMACS .top file detected by MDAnalysis. Passing ParmEd object instead...")
        try:
            u = mda.Universe(top, traj_file)
            ref_u = mda.Universe(top, traj_file)
        except UnicodeDecodeError:
            u = mda.Universe(top, traj_file, format='NCDF')
            ref_u = mda.Universe(top, traj_file, format='NCDF')
            forced_ncdf = True
    except UnicodeDecodeError:
        print("  -> UnicodeDecodeError detected. Assuming binary NetCDF format. Forcing NCDF reader...")
        u = mda.Universe(top_file, traj_file, format='NCDF')
        ref_u = mda.Universe(top_file, traj_file, format='NCDF')
        forced_ncdf = True

    t_start = params['traj'][0] - 1
    t_step = params['traj'][2]
    t_stop = None if params['traj'][1] == 'last' else params['traj'][1]
    
    ref_u.trajectory[t_start]
    
    fit_opt = params['rmsfit']
    fit_sel = "none"
    if fit_opt == 'host':    
        fit_sel = get_mda_selection(params['host'])
    elif fit_opt == 'guest': 
        fit_sel = get_mda_selection(params['guest'])
    elif fit_opt == 'complex': 
        fit_sel = get_mda_selection(params['host'] + params['guest'])

    if params['centering']:
        if fit_sel != "none":
            ref_com = ref_u.select_atoms(fit_sel).center_of_mass()
        else:
            ref_com = ref_u.select_atoms(get_mda_selection(params['host'])).center_of_mass()
        ref_u.atoms.translate(-ref_com)

    writer = None
    new_traj_file = None
    if params['newtraj']:
        base, ext = os.path.splitext(traj_file)
        new_traj_file = f"{base}.new{ext}"
        print(f"Setting up new trajectory output: {new_traj_file}")
        
        try:
            if forced_ncdf:
                writer = mda.Writer(new_traj_file, u.atoms.n_atoms, format='NCDF')
            else:
                writer = mda.Writer(new_traj_file, u.atoms.n_atoms)
        except Exception as e:
            print(f"Warning: Failed to initialize trajectory writer. Error: {e}")
            writer = None

    # --- Step 4: Iterating over Trajectory Frames ---
    results = []
    
    print("\nStarting trajectory analysis loop...")
    for ts in u.trajectory[t_start:t_stop:t_step]:
        current_frame = ts.frame + 1 
        print(f"\n--- Processing Frame {current_frame} ---")
        
        if fit_sel != "none":
            align.alignto(u, ref_u, select=fit_sel, match_atoms=True)
        elif params['centering']:
            current_com = u.select_atoms(get_mda_selection(params['host'])).center_of_mass()
            u.atoms.translate(-current_com)
            
        if writer is not None:
            writer.write(u.atoms)
            
        # E_MM is calculated using original relative positions BEFORE individual state centering
        if len(guest_idx_0) > 0:
            E_MM = calc_emm(
                u.atoms[host_idx_0].positions, u.atoms[guest_idx_0].positions,
                sigmas[host_idx_0], sigmas[guest_idx_0],
                epsilons[host_idx_0], epsilons[guest_idx_0],
                charges[host_idx_0], charges[guest_idx_0]
            )
        else:
            E_MM = 0.0
            
        components = {'Frame': current_frame, 'E_MM': E_MM}
        
        c_flag = params['centereach']
        
        # --- Process Host State ---
        h_inp = f"host_{current_frame}.inp"
        h_out = f"host_{current_frame}.out"
        h_xmu = f"host_{current_frame}.xmu"
        
        print("  Running RISMiCal for Host...")
        write_rism_inp(h_inp, base_inp, u, host_idx_0, names, sigmas, epsilons, charges, f"host_{current_frame}", centereach=c_flag)
        subprocess.run(f"rismical.x 3d {h_inp} > {h_out}", shell=True, check=True)
        components['sfe_h'], components['corr_h'] = extract_xmu(h_xmu)
        
        # --- Process Guest State ---
        if len(guest_idx_0) > 0:
            g_inp = f"guest_{current_frame}.inp"
            g_out = f"guest_{current_frame}.out"
            g_xmu = f"guest_{current_frame}.xmu"
            
            print("  Running RISMiCal for Guest...")
            write_rism_inp(g_inp, base_inp, u, guest_idx_0, names, sigmas, epsilons, charges, f"guest_{current_frame}", centereach=c_flag)
            subprocess.run(f"rismical.x 3d {g_inp} > {g_out}", shell=True, check=True)
            components['sfe_g'], components['corr_g'] = extract_xmu(g_xmu)
        else:
            components['sfe_g'], components['corr_g'] = 0.0, 0.0
            
        # --- Process Complex State ---
        if len(guest_idx_0) > 0:
            c_inp = f"complex_{current_frame}.inp"
            c_out = f"complex_{current_frame}.out"
            c_xmu = f"complex_{current_frame}.xmu"
            
            print("  Running RISMiCal for Complex...")
            write_rism_inp(c_inp, base_inp, u, complex_idx_0, names, sigmas, epsilons, charges, f"complex_{current_frame}", centereach=c_flag)
            subprocess.run(f"rismical.x 3d {c_inp} > {c_out}", shell=True, check=True)
            components['sfe_c'], components['corr_c'] = extract_xmu(c_xmu)
        else:
            components['sfe_c'], components['corr_c'] = components['sfe_h'], components['corr_h']

        components['dSFE'] = components['sfe_c'] - components['sfe_h'] - components['sfe_g']
        components['dPC'] = components['corr_c'] - components['corr_h'] - components['corr_g']
        components['dG_bind'] = components['E_MM'] + components['dSFE'] + components['dPC']
        
        results.append(components)
        print(f"  Frame {current_frame} dG_bind = {components['dG_bind']/1000:.4f} kJ/mol")

    if writer is not None:
        writer.close()
        print(f"\nNew aligned trajectory has been saved to: {new_traj_file}")

    if results:
        base_inp_name, _ = os.path.splitext(inp_file)
        csv_filename = f"{base_inp_name}.csv"
        print(f"Exporting frame-by-frame data to: {csv_filename}")
        
        csv_headers = [
            'Frame', 'E_MM', 
            'sfe_c', 'sfe_h', 'sfe_g', 'dSFE', 
            'corr_c', 'corr_h', 'corr_g', 'dPC', 
            'dG_bind'
        ]
        
        with open(csv_filename, mode='w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(csv_headers)
            for row in results:
                csv_writer.writerow([row[key] for key in csv_headers])

    if not results:
        print("No frames were processed. Exiting.")
        return

    keys = ['E_MM', 'dSFE', 'dPC', 'dG_bind']
    stats = {k: [r[k] for r in results] for k in keys}
    
    print("\n" + "="*50)
    print("      Final Free Energy Summary (kJ/mol)")
    print("="*50)
    print(f"{'Component':<12} | {'Mean':>15} | {'Variance':>15}")
    print("-" * 50)
    
    for k in keys:
        mean_val = np.mean(stats[k])
        var_val = np.var(stats[k], ddof=1) if len(stats[k]) > 1 else 0.0
        print(f"{k:<12} | {mean_val/1000:15.4f} | {var_val/1000000:15.4f}")
    print("="*50)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python rismical-bind.py <topology_file> <trajectory_file> <base_inp_file>")
        sys.exit(1)
        
    main(sys.argv[1], sys.argv[2], sys.argv[3])
