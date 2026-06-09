#!/usr/bin/env python3
"""
RISMiCal-QM-PySCF.py
An open-source, purely in-memory wrapper for QM/3D-RISM calculations using PySCF.
Includes support for Picard and DIIS methods, a highly memory-optimized 
Electrostatic Potential (MEP) generator for large 3D-RISM grids, Franck-Condon 
state calculations, and automatic Molden file generation for full MO visualization.
"""

import os
import sys
import shutil
import subprocess
import re
import numpy as np
from scipy.spatial.distance import cdist

try:
    from pyscf import gto, dft, qmmm, tdscf, lib, df
    from pyscf.tools import cubegen, molden
except ImportError:
    print("[ERROR] PySCF is not installed. Please run: pip install pyscf")
    sys.exit(1)

# ==============================================================================
# Commands and Paths Configuration
# ==============================================================================
RISMICAL_CMD = "rismical.x"

# Physical constants for energy and length conversions
HARTREE_TO_JMOL = 2625499.6394799
COULOMB_TO_JMOL = 1389354.56      # 1 e^2 / Angstrom -> J/mol
ANG_TO_BOHR     = 1.8897261246    # Angstrom -> Bohr

def parse_fortran_float(val_str):
    try: return float(str(val_str).lower().replace('d', 'e'))
    except ValueError: return 0.0

def parse_qmpart(qmpart_str, total_atoms):
    if not qmpart_str or not str(qmpart_str).strip():
        return list(range(1, total_atoms + 1))
    indices = []
    for p in str(qmpart_str).split(','):
        if not p.strip(): continue
        val = int(p.strip())
        if val > 0:
            indices.append(val)
            last_idx = val
        elif val < 0:
            indices.extend(range(last_idx + 1, abs(val) + 1))
            last_idx = abs(val)
    return sorted(list(set(indices)))

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
        for j, p in enumerate(placeholders): val_str = val_str.replace(f"__QUOTE_{j}__", p)
        if (val_str.startswith('"') and val_str.endswith('"')) or (val_str.startswith("'") and val_str.endswith("'")):
            val_str = val_str[1:-1]
        params[key] = val_str

def read_input_file(inp_file):
    params, udata, udata_line_indices = {}, [], []
    with open(inp_file, 'r', errors='replace') as f: lines = f.readlines()
    in_rismicalqm, in_grid3d, in_udata = False, False, False
    rismicalqm_block_str, grid3d_block_str = "", ""
    for i, line in enumerate(lines):
        clean_line = line.split('!')[0].strip()
        if not clean_line: continue
        upper_line = clean_line.upper()
        if upper_line.startswith('$RISMICALQM') or upper_line.startswith('&RISMICALQM'):
            in_rismicalqm = True
            header_len = len('$RISMICALQM') if upper_line.startswith('$RISMICALQM') else len('&RISMICALQM')
            rismicalqm_block_str += clean_line[header_len:] + " "
        elif upper_line.startswith('$GRID3D') or upper_line.startswith('&GRID3D'):
            in_grid3d = True
            header_len = len('$GRID3D') if upper_line.startswith('$GRID3D') else len('&GRID3D')
            grid3d_block_str += clean_line[header_len:] + " "
        elif upper_line.startswith('$UDATA') or upper_line.startswith('&UDATA'):
            in_udata = True
        elif in_rismicalqm and (upper_line == '$END' or upper_line == '/'): in_rismicalqm = False
        elif in_grid3d and (upper_line == '$END' or upper_line == '/'): in_grid3d = False
        elif in_udata and (upper_line == '$END' or upper_line == '/'): in_udata = False
        elif in_rismicalqm: rismicalqm_block_str += clean_line + " "
        elif in_grid3d: grid3d_block_str += clean_line + " "
        elif in_udata:
            parts = clean_line.split()
            if len(parts) >= 7: udata.append(parts); udata_line_indices.append(i)
    extract_namelist(rismicalqm_block_str, params)
    extract_namelist(grid3d_block_str, params)
    return params, udata, udata_line_indices, lines

def build_pyscf_mol(params, udata, qm_indices):
    mol_str = ""
    for idx in qm_indices:
        atom = udata[idx-1]
        x, y, z = parse_fortran_float(atom[4]), parse_fortran_float(atom[5]), parse_fortran_float(atom[6])
        mol_str += f"{atom[0]} {x} {y} {z}; "
    
    basis  = params.get('basis', '6-31g*')
    charge = int(params.get('charge', '0'))
    spin   = int(params.get('spin', '0'))
    
    mol = gto.Mole()
    mol.atom = mol_str
    mol.basis = basis
    mol.charge = charge
    mol.spin = spin
    mol.verbose = 0
    mol.build()
    return mol

# [MODIFIED] Added mm_coords and mm_charges explicitly to separate E_UU and E_QM-MM properly
def run_pyscf_scf(mol, params, ext_coords=None, ext_charges=None, mm_coords=None, mm_charges=None):
    xc = params.get('xc', 'b3lyp')
    mf = dft.RKS(mol)
    mf.xc = xc
    
    if ext_coords is not None and len(ext_coords) > 0:
        # [FIXED] Explicitly set unit='Angstrom'. Without this, PySCF defaults to Bohr, 
        # causing a massive 1.889x error in external electrostatic interactions.
        mf = qmmm.mm_charge(mf, ext_coords, ext_charges, unit='Angstrom')
        
    mf.kernel()
    if not mf.converged:
        print("\n  [ERROR] PySCF SCF calculation did not converge!")
        sys.exit(1)
        
    dm = mf.make_rdm1()
    
    # --- Energy Breakdown to match G16 definition (E_UU, E_QM-MM) ---
    mf_base = dft.RKS(mol)
    mf_base.xc = xc
    # E_UU: Pure polarized QM internal energy without external field
    e_uu_elec = mf_base.energy_elec(dm)[0]
    e_uu_h = e_uu_elec + mol.energy_nuc()
    
    # E_QM-MM: Interaction strictly between QM density/nuclei and MM atoms
    e_qm_mm_h = 0.0
    if mm_coords is not None and len(mm_coords) > 0:
        mf_mm = dft.RKS(mol)
        mf_mm = qmmm.mm_charge(mf_mm, mm_coords, mm_charges, unit='Angstrom')
        h_ext_mm = mf_mm.get_hcore() - mf_base.get_hcore()
        e_qm_mm_elec = np.einsum('ij,ji->', h_ext_mm, dm)
        
        e_qm_mm_nuc = 0.0
        mm_c_bohr = mm_coords * ANG_TO_BOHR 
        for i in range(mol.natm):
            q = mol.atom_charge(i)
            r = mol.atom_coord(i) # PySCF coordinates are in Bohr
            dists = np.linalg.norm(mm_c_bohr - r, axis=1)
            dists = np.where(dists < 1e-12, 1e-12, dists)
            e_qm_mm_nuc += np.sum(q * mm_charges / dists)
            
        e_qm_mm_h = e_qm_mm_elec + e_qm_mm_nuc
    
    # ----------------------------------------------------------------
    
    td_energies_h = []
    if 'td' in params or 'cis' in params:
        td_obj = tdscf.TDDFT(mf)
        td_obj.nstates = int(params.get('nstates', '3'))
        td_obj.kernel()
        target_root = int(params.get('root', '1'))
        if target_root > len(td_obj.e):
            print(f"\n  [ERROR] Requested root {target_root} exceeds computed states.")
            sys.exit(1)
            
        td_energies_h = td_obj.e 
        
        X, Y = td_obj.xy[target_root - 1]
        mo_coeff = mf.mo_coeff
        mo_occ = mf.mo_occ
        nocc = np.count_nonzero(mo_occ > 0)
        
        dP_mo = np.zeros((mo_occ.size, mo_occ.size))
        dP_mo[:nocc, :nocc] = -2.0 * (X @ X.T + Y @ Y.T)
        dP_mo[nocc:, nocc:] =  2.0 * (X.T @ X + Y.T @ Y)
        
        dm_ex_ao = mo_coeff @ dP_mo @ mo_coeff.T
        dm = mf.make_rdm1() + dm_ex_ao
        
        _, qm_charges = mf.mulliken_pop(mol=mol, dm=dm, s=mf.get_ovlp(), verbose=0)
        
        # Re-evaluate E_UU and E_QM-MM for the excited state density
        e_uu_elec = mf_base.energy_elec(dm)[0]
        e_uu_h = e_uu_elec + mol.energy_nuc()
        if mm_coords is not None and len(mm_coords) > 0:
            e_qm_mm_elec = np.einsum('ij,ji->', h_ext_mm, dm)
            e_qm_mm_h = e_qm_mm_elec + e_qm_mm_nuc
    else:
        _, qm_charges = mf.mulliken_pop(verbose=0)
        
    return e_uu_h, e_qm_mm_h, td_energies_h, dm, qm_charges, mf

def process_qv_and_get_ext_charges(qv_file, udata, qvcutoff, qvcore, mm_indices):
    ext_c, ext_q = [], []
    e_mv_jmol = 0.0
    if os.path.exists(qv_file):
        qv_coords, qv_charges = [], []
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
                        
        if qv_coords:
            qv_c_arr, qv_q_arr = np.array(qv_coords), np.array(qv_charges)
            all_solute_coords = [[parse_fortran_float(row[4]), parse_fortran_float(row[5]), parse_fortran_float(row[6])] for row in udata]
            
            dists = cdist(qv_c_arr, np.array(all_solute_coords))
            valid_mask = np.min(dists, axis=1) >= qvcore
            solv_c, solv_q = qv_c_arr[valid_mask], qv_q_arr[valid_mask]
            
            ext_c.extend(solv_c.tolist())
            ext_q.extend(solv_q.tolist())

            if mm_indices and len(solv_c) > 0:
                mm_c_arr = np.array([[parse_fortran_float(udata[i-1][4]), parse_fortran_float(udata[i-1][5]), parse_fortran_float(udata[i-1][6])] for i in mm_indices])
                mm_q_arr = np.array([parse_fortran_float(udata[i-1][3]) for i in mm_indices])
                inv_dists = np.where(cdist(mm_c_arr, solv_c) > 1e-6, 1.0/cdist(mm_c_arr, solv_c), 0.0)
                e_mv_jmol = np.dot(mm_q_arr, np.dot(inv_dists, solv_q)) * COULOMB_TO_JMOL

    if mm_indices:
        for i in mm_indices:
            atom = udata[i-1]
            ext_c.append([parse_fortran_float(atom[4]), parse_fortran_float(atom[5]), parse_fortran_float(atom[6])])
            ext_q.append(parse_fortran_float(atom[3]))
            
    if ext_c:
        with open(qv_file, 'w') as fout:
            for c, q in zip(ext_c, ext_q):
                fout.write(f" {c[0]:11.6f} {c[1]:11.6f} {c[2]:11.6f}  {q:.6e}\n")

    return np.array(ext_c) if ext_c else None, np.array(ext_q) if ext_q else None, e_mv_jmol

def generate_and_write_cube(mol, dm, outfile, ngrid3d, rdelta3d, udata, mm_indices, mix_method='picard', mix_alpha=1.0, prev_pot=None, diis_engine=None):
    """
    Computes MEP directly in memory in a chunked manner, adds MM potential, 
    applies mixing (Picard or DIIS) to stabilize convergence, 
    and writes the final 3D-RISM ready Cube file. Returns the mixed potential array.
    """
    origin_ang = -rdelta3d * (ngrid3d / 2.0)
    origin_bohr = origin_ang * ANG_TO_BOHR
    step_bohr = rdelta3d * ANG_TO_BOHR
    extent_bohr = step_bohr * (ngrid3d - 1)
    
    cube = cubegen.Cube(mol, nx=ngrid3d, ny=ngrid3d, nz=ngrid3d)
    cube.boxorig = np.array([origin_bohr, origin_bohr, origin_bohr])
    cube.box = np.diag([extent_bohr, extent_bohr, extent_bohr])
    
    coords = cube.get_coords()
    npts = coords.shape[0]
    
    if dm.ndim == 3:
        dm_tot = dm[0] + dm[1]
    else:
        dm_tot = dm
        
    qm_pot_flat = np.zeros(npts)
    
    for i in range(mol.natm):
        q = mol.atom_charge(i)
        r = mol.atom_coord(i)
        d = np.linalg.norm(coords - r, axis=1)
        qm_pot_flat += q / np.maximum(d, 1e-12)
        
    nbas = mol.nao_nr()
    chunk_size = int((500 * 1024**2) / (nbas * nbas * 8)) 
    chunk_size = max(1000, min(chunk_size, npts))
    
    for p0 in range(0, npts, chunk_size):
        p1 = min(p0 + chunk_size, npts)
        fakemol = gto.fakemol_for_charges(coords[p0:p1])
        int3c = df.incore.aux_e2(mol, fakemol, intor='int3c2e', aosym='s1')
        v_elec = np.einsum('ij,ijx->x', dm_tot, int3c)
        qm_pot_flat[p0:p1] -= v_elec
        
    pot = qm_pot_flat.reshape(ngrid3d, ngrid3d, ngrid3d)
    
    if mm_indices:
        X = cube.boxorig[0] + np.arange(ngrid3d) * step_bohr
        Y = cube.boxorig[1] + np.arange(ngrid3d) * step_bohr
        Z = cube.boxorig[2] + np.arange(ngrid3d) * step_bohr
        XX, YY, ZZ = np.meshgrid(X, Y, Z, indexing='ij')
        XX_flat, YY_flat, ZZ_flat = XX.ravel(), YY.ravel(), ZZ.ravel()
        
        v_mm = np.zeros_like(XX_flat)
        for i in mm_indices:
            atom = udata[i-1]
            c_bohr = np.array([parse_fortran_float(atom[4]), parse_fortran_float(atom[5]), parse_fortran_float(atom[6])]) * ANG_TO_BOHR
            q = parse_fortran_float(atom[3])
            
            d_sq = (XX_flat - c_bohr[0])**2 + (YY_flat - c_bohr[1])**2 + (ZZ_flat - c_bohr[2])**2
            v_mm += q / np.sqrt(np.maximum(d_sq, 1e-12))
            
        pot += v_mm.reshape(ngrid3d, ngrid3d, ngrid3d)
        
    new_pot = pot
    mixed_pot = new_pot

    if prev_pot is not None:
        if mix_method == 'diis' and diis_engine is not None:
            print("  [INFO] Applying DIIS extrapolation for electrostatic potential.")
            err_vec = new_pot - prev_pot
            try:
                mixed_pot = diis_engine.update(new_pot, xerr=err_vec)
            except Exception as e:
                print(f"  [WARNING] DIIS update failed: {e}. Falling back to Picard mixing (alpha={mix_alpha}).")
                mixed_pot = mix_alpha * new_pot + (1.0 - mix_alpha) * prev_pot
                
        elif mix_method == 'picard' and mix_alpha < 1.0:
            print(f"  [INFO] Applying Picard mixing with alpha = {mix_alpha}")
            mixed_pot = mix_alpha * new_pot + (1.0 - mix_alpha) * prev_pot
            
    cube.write(mixed_pot, outfile, comment="QM+MM Electrostatic Potential generated by RISMiCal-PySCF")
    return mixed_pot

def generate_molden_file(mf, base_name):
    """
    Writes all Molecular Orbitals, energies, and occupancies to a single Molden file 
    for comprehensive visualization in tools like Avogadro.
    """
    molden_file = f"{base_name}.molden"
    print(f"\n--- Generating Molden File for Avogadro ---")
    print(f"  Writing all MOs to {molden_file}...")
    molden.from_scf(mf, molden_file)

def update_rismical_input(inp_file, lines, udata, udata_line_indices, qm_indices, charges):
    for c_idx, idx in enumerate(qm_indices):
        row = udata[idx-1]
        p1, p2 = parse_fortran_float(row[1]), parse_fortran_float(row[2])
        x, y, z = parse_fortran_float(row[4]), parse_fortran_float(row[5]), parse_fortran_float(row[6])
        lines[udata_line_indices[idx-1]] = f"{row[0]:<7s}{p1:8.4f}{p2:10.4f}{charges[c_idx]:12.6f}     {x:8.4f}   {y:8.4f}   {z:8.4f}\n"
    with open(inp_file, 'w', errors='replace') as f: f.writelines(lines)

def check_rismical_termination(rsmout_file):
    if not os.path.exists(rsmout_file): return False
    with open(rsmout_file, 'r', errors='replace') as f:
        for line in f.readlines()[-30:]:
            if "RISMiCal computation is completed normally" in line: return True
    return False

def read_xmu(xmu_file):
    sfe_sc, se_es = 0.0, 0.0
    with open(xmu_file, 'r', errors='replace') as f:
        for line in f:
            if 'SFE_SC=' in line: sfe_sc = float(line.split('=')[1].split('!')[0].strip())
            elif 'SE_ES=' in line: se_es = float(line.split('=')[1].split('!')[0].strip())
    return sfe_sc, se_es

def main():
    fc_mode = False
    if "-FC" in sys.argv:
        fc_mode = True
        sys.argv.remove("-FC")

    if len(sys.argv) < 2: 
        print("Usage: python RISMiCal-QM-PySCF.py <input_file> [-FC]"); sys.exit(1)
        
    inp_file = sys.argv[1]
    base_name = os.path.splitext(inp_file)[0]
    ascii_cube = f"{base_name}.cub"
    xmu_file, qv_file, org_inp_file = f"{base_name}.xmu", f"{base_name}.qv", f"{base_name}.org.inp"
    
    if not os.path.exists(org_inp_file): shutil.copy(inp_file, org_inp_file)
    
    params, udata, udata_indices, all_lines = read_input_file(inp_file)
    
    xc = params.get('xc', 'b3lyp')
    basis = params.get('basis', '6-31g*')
    
    scfconv  = parse_fortran_float(params.get('scfconv', '1e-4'))
    qvcutoff = parse_fortran_float(params.get('qvcutoff', '1e-6'))
    qvcore   = parse_fortran_float(params.get('qvcore', '0.5'))
    ngrid3d  = int(parse_fortran_float(params.get('ngrid3d', '128')))
    rdelta3d = parse_fortran_float(params.get('rdelta3d', '0.5'))
    
    mix_method = params.get('mix_method', 'picard').lower()
    mix_alpha  = parse_fortran_float(params.get('mix_alpha', '1.0'))
    
    diis_engine = None
    if mix_method == 'diis':
        diis_engine = lib.diis.DIIS()
        diis_engine.space = 8
    
    total_atoms = len(udata)
    qm_indices = parse_qmpart(params.get('qmpart', ''), total_atoms)
    mm_indices = [i for i in range(1, total_atoms + 1) if i not in qm_indices]
    
    # Pre-calculate explicit arrays for MM to accurately compute E_QM-MM later
    mm_c, mm_q = None, None
    if mm_indices:
        mm_c = np.array([[parse_fortran_float(udata[i-1][4]), parse_fortran_float(udata[i-1][5]), parse_fortran_float(udata[i-1][6])] for i in mm_indices])
        mm_q = np.array([parse_fortran_float(udata[i-1][3]) for i in mm_indices])
    
    print(f"--- 3D-RISM-PySCF Started ({base_name}) ---")
    print(f" Engine     : PySCF ({xc}/{basis})")
    print(f" Total Atoms: {total_atoms} (QM: {len(qm_indices)}, MM: {len(mm_indices)})")
    
    if mix_method == 'diis':
        print(" Convergence Acceleration: DIIS Method")
    elif mix_alpha < 1.0:
        print(f" Convergence Acceleration: Picard Mixing (alpha = {mix_alpha})")
    
    mol = build_pyscf_mol(params, udata, qm_indices)

    if fc_mode:
        print("\n--- Franck-Condon (FC) State Calculation ---")
        if not os.path.exists(qv_file):
            print(f"  [ERROR] FC mode requires an existing '{qv_file}'.")
            sys.exit(1)
            
        print(f"  Using existing frozen solvent/MM charges from: {qv_file}")
        ext_coords, ext_charges, e_mv = process_qv_and_get_ext_charges(qv_file, udata, qvcutoff, qvcore, mm_indices)
        
        e_uu_h, e_qm_mm_h, td_energies_h, dm, new_chg, final_mf = run_pyscf_scf(mol, params, ext_coords, ext_charges, mm_c, mm_q)
        
        e_uu = e_uu_h * HARTREE_TO_JMOL / 1000.0
        e_qm_mm = e_qm_mm_h * HARTREE_TO_JMOL / 1000.0
        e_mv_kj = e_mv / 1000.0
        
        generate_molden_file(final_mf, base_name)
        
        print("\n=======================================================")
        print(" SUMMARY OF RISMiCal-QM (FC MODE) CALCULATION")
        print("=======================================================")
        print(" >>> FC Calculation Completed Successfully! <<<")
        print(f" QM Internal Energy (E_UU)   : {e_uu:15.5f} kJ/mol")
        if mm_indices:
            print(f" QM-MM Int. (E_QM-MM)        : {e_qm_mm:15.5f} kJ/mol")
            
        if np.asarray(td_energies_h).size > 0:
            print("-------------------------------------------------------")
            print(" Excited States Summary:")
            target_root = int(params.get('root', '1'))
            for i, e_h in enumerate(td_energies_h):
                e_ev = e_h * 27.211386245988 
                marker = " (*Target)" if (i + 1) == target_root else ""
                print(f"   State {i+1:2d} : {e_ev:15.5f} eV{marker}")
                
        print("-------------------------------------------------------")
        print(" Effective Charges (ESP/MK) for QM region:")
        for i, q in zip(qm_indices, new_chg):
            atom_name = udata[i-1][0]
            print(f"   Atom {i:2d} ({atom_name:2s}) : {q:9.6f} e")
        print("-------------------------------------------------------")
        print(" Generated Files:")
        print(f"   Input Backup  : {org_inp_file}")
        print(f"   Molden File   : {base_name}.molden")
        print(f"   Ext. Charges  : {qv_file}")
        print("=======================================================")
        sys.exit(0) 

    # ---------------------------------------------------------
    # Pre-Step 1: Pure QM Vacuum
    # ---------------------------------------------------------
    print("\n--- Pre-Step 1: QM Vacuum (E_gas) ---")
    e_uu_h, _, _, _, _, _ = run_pyscf_scf(mol, params)
    e_gas = e_uu_h * HARTREE_TO_JMOL / 1000.0
    
    # ---------------------------------------------------------
    # Pre-Step 2: QM + MM Vacuum
    # ---------------------------------------------------------
    print("--- Pre-Step 2: QM+MM Vacuum (E_QMMM_gas) ---")
    if mm_indices:
        e_uu_h, e_qm_mm_h, _, _, _, _ = run_pyscf_scf(mol, params, ext_coords=mm_c, ext_charges=mm_q, mm_coords=mm_c, mm_charges=mm_q)
        e_qmmm_gas = (e_uu_h + e_qm_mm_h) * HARTREE_TO_JMOL / 1000.0
    else:
        e_qmmm_gas = e_gas
        print("  No MM atoms specified. Skipping...")
        
    e_qmmm_int = e_qmmm_gas - e_gas

    # ---------------------------------------------------------
    # Initializing Solvent Distribution
    # ---------------------------------------------------------
    print("\n--- Initializing Solvent Distribution ---")
    _, _, _, dm_qmmm, init_qm_charges, _ = run_pyscf_scf(mol, params, ext_coords=mm_c, ext_charges=mm_q, mm_coords=mm_c, mm_charges=mm_q)
    
    prev_pot = generate_and_write_cube(mol, dm_qmmm, ascii_cube, ngrid3d, rdelta3d, udata, mm_indices, 
                                       mix_method='none', mix_alpha=1.0, prev_pot=None, diis_engine=None)
                                       
    update_rismical_input(inp_file, all_lines, udata, udata_indices, qm_indices, init_qm_charges)
    current_rsmout = f"{base_name}.rsmout.0"
    with open(current_rsmout, "w", errors='replace') as f_rsm: 
        subprocess.run([RISMICAL_CMD, "3d", inp_file], stdout=f_rsm, stderr=subprocess.STDOUT)
    if not check_rismical_termination(current_rsmout):
        print(f"  [ERROR] RISMiCal did not terminate normally. Check {current_rsmout}"); sys.exit(1)

    # ---------------------------------------------------------
    # Main SCF Loop
    # ---------------------------------------------------------
    prev_g = 0.0
    iter_num = 1
    final_mf = None
    td_energies_h = []
    print("\n--- Starting QM/MM/3D-RISM SCF Loop ---")
    while True:
        print(f"\n[Iteration {iter_num}]")
        
        sfe, se = read_xmu(xmu_file)
        ext_coords, ext_charges, e_mv = process_qv_and_get_ext_charges(qv_file, udata, qvcutoff, qvcore, mm_indices)
        
        e_uu_h, e_qm_mm_h, td_energies_h, dm, new_chg, final_mf = run_pyscf_scf(mol, params, ext_coords, ext_charges, mm_c, mm_q)
        
        e_uu = e_uu_h * HARTREE_TO_JMOL / 1000.0
        e_qm_mm = e_qm_mm_h * HARTREE_TO_JMOL / 1000.0
        se_kj = se / 1000.0
        sfe_kj = sfe / 1000.0
        
        g_tot = e_uu + e_qm_mm + sfe_kj
        
        print(f"  E_UU    = {e_uu:.4f} kJ/mol")
        if mm_indices: print(f"  E_QM-MM = {e_qm_mm:.4f} kJ/mol")
        print(f"  E_QMMM-v= {se_kj:.4f} kJ/mol")
        print(f"  SFE_SC  = {sfe_kj:.4f} kJ/mol")
        print(f"  G_TOTAL = {g_tot:.4f} kJ/mol")
        
        d_e_h = abs(g_tot - prev_g) / (HARTREE_TO_JMOL / 1000.0)
        print(f"  Delta G = {d_e_h:.6e} Hartree (Threshold: {scfconv})")
        
        if iter_num > 1 and d_e_h <= scfconv: break
        prev_g = g_tot
        
        prev_pot = generate_and_write_cube(mol, dm, ascii_cube, ngrid3d, rdelta3d, udata, mm_indices,
                                           mix_method=mix_method, mix_alpha=mix_alpha, prev_pot=prev_pot, diis_engine=diis_engine)
                                           
        update_rismical_input(inp_file, all_lines, udata, udata_indices, qm_indices, new_chg)
        current_rsmout = f"{base_name}.rsmout.{iter_num}"
        with open(current_rsmout, "w", errors='replace') as f_rsm: 
            subprocess.run([RISMICAL_CMD, "3d", inp_file], stdout=f_rsm, stderr=subprocess.STDOUT)
            
        if not check_rismical_termination(current_rsmout):
            print(f"  [ERROR] RISMiCal did not terminate normally in Iteration {iter_num}. Check {current_rsmout}"); sys.exit(1)
            
        iter_num += 1
        
    if final_mf is not None:
        generate_molden_file(final_mf, base_name)

    print("\n=======================================================")
    print(" SUMMARY OF RISMiCal-QM CALCULATION")
    print("=======================================================")
    print(" >>> SCF Converged! <<<")
    print(f" QM Gas Energy (E_QM_gas)    : {e_gas:15.5f} kJ/mol")
    if mm_indices:
        print(f" QMMM Gas Energy (E_QMMM_gas): {e_qmmm_gas:15.5f} kJ/mol")
        print(f" QM-MM Interaction (gas)     : {e_qmmm_int:15.5f} kJ/mol")
    print("-" * 55)
    print(f" Total Free Energy (G_TOTAL) : {g_tot:15.5f} kJ/mol")
    print(f" QM Internal Energy (E_UU)   : {e_uu:15.5f} kJ/mol")
    if mm_indices:
        print(f" QM-MM Int. (E_QM-MM)        : {e_qm_mm:15.5f} kJ/mol")
    print(f" QM/MM-Solvent Int. (E_QMMM-v): {se_kj:15.5f} kJ/mol")
    print(f" Solvation Free E. (SFE_SC)  : {sfe_kj:15.5f} kJ/mol")
    print("-" * 55)
    
    if np.asarray(td_energies_h).size > 0:
        print(" Excited States Summary:")
        target_root = int(params.get('root', '1'))
        for i, e_h in enumerate(td_energies_h):
            e_ev = e_h * 27.211386245988 
            marker = " (*Target)" if (i + 1) == target_root else ""
            print(f"   State {i+1:2d} : {e_ev:15.5f} eV{marker}")
        print("-" * 55)
        
    print(" Effective Charges (ESP/MK) for QM region:")
    for i, q in zip(qm_indices, new_chg):
        atom_name = udata[i-1][0]
        print(f"   Atom {i:2d} ({atom_name:2s}) : {q:9.6f} e")
    print("-" * 55)
    print(" Generated Files:")
    print(f"   Input Backup  : {org_inp_file}")
    print(f"   RISMiCal LOGs : {base_name}.rsmout.*")
    print(f"   ASCII Cube    : {ascii_cube}")
    print(f"   Molden File   : {base_name}.molden")
    print(f"   RISMiCal OUT  : {xmu_file}")
    print(f"   Ext. Charges  : {qv_file}")
    print("=======================================================")

if __name__ == "__main__": 
    main()