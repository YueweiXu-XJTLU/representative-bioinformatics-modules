############################################
# 0. REINITIALIZE & LOAD NATIVE & DOCKED STRUCTURES
############################################
reinitialize everything

select rec_native, native and chain A and polymer.protein
select lig_native, native and chain E and polymer.protein

load ../cluspro_results/model.000.00.pdb
set_name rec.pdb, rec_dock_raw
set_name lig.000.00.pdb, lig_dock_raw
create dock, rec_dock_raw or lig_dock_raw
create ace2_rbd, dock and chain A+E and polymer.protein

select ace2, ace2_rbd and chain A
select rbd,  ace2_rbd and chain E

hide everything
show cartoon, ace2_rbd
color lightblue, ace2
color tv_orange, rbd

############################################################
# METHOD 1: Pymol
############################################################

# 1. CONTACT-BASED (≤ 4 Å)
select ace2_contact, byres (ace2 within 4.0 of rbd)
select rbd_contact,  byres (rbd  within 4.0 of ace2)

create 1_1_res, ace2_contact or rbd_contact
dist   1_1_bonds, ace2_contact, rbd_contact, 4.0
show   sticks, 1_1_res
show   dash,   1_1_bonds
color yellow, 1_1_res
color yellow, 1_1_bonds

# 2. BSA-BASED (ΔSASA) — PyMOL approximation of PISA
# Make monomers
create ace2_monomer, ace2
create rbd_monomer,  rbd

# Write SASA into b-factor
get_area ace2_monomer, load_b=1
get_area rbd_monomer,  load_b=1
get_area ace2, load_b=1
get_area rbd,  load_b=1

# Compute ΔSASA for each residue
python
from pymol import cmd

def get_deltasasa(free_sel, complex_sel):
    free = cmd.get_model(free_sel).atom
    comp = cmd.get_model(complex_sel).atom
    ds = {}
    for a, b in zip(free, comp):
        key = (a.resi, a.resn)
        ds.setdefault(key, 0)
        ds[key] += a.b - b.b  # ΔSASA = SASA_monomer − SASA_complex
    return ds

ds_ace2 = get_deltasasa("ace2_monomer", "ace2")
ds_rbd  = get_deltasasa("rbd_monomer",  "rbd")

# Identify interface residues as ΔSASA > 1 Å²
ace2_bsa_res = [resi for (resi,resn),v in ds_ace2.items() if v > 1]
rbd_bsa_res  = [resi for (resi,resn),v in ds_rbd.items()  if v > 1]

print("ACE2 ΔSASA residues:", ace2_bsa_res)
print("RBD  ΔSASA residues:", rbd_bsa_res)

# Create PyMOL selections
cmd.select("ace2_bsa", "ace2 and resi " + "+".join(ace2_bsa_res))
cmd.select("rbd_bsa",  "rbd and resi "  + "+".join(rbd_bsa_res))
python end

create 1_2_res, ace2_bsa or rbd_bsa
show   sticks, 1_2_res
color green, 1_2_res

# 3. HYDROGEN BONDS / POLAR CONTACTS               #
# Use PyMOL hydrogen-bond detector (mode=2)
dist hbonds, ace2, rbd, mode=2

select ace2_hbond, byres (ace2 within 3.6 of rbd)
select rbd_hbond,  byres (rbd  within 3.6 of ace2)

create 1_3_res, ace2_hbond or rbd_hbond
dist   1_3_bonds, ace2_hbond, rbd_hbond, 3.6
show   sticks, 1_3_res
show   dash,   1_3_bonds
color red, 1_3_res
color red, 1_3_bonds

# 4. ELECTROSTATIC SALT BRIDGES (≤ 4 Å)
# Positive ACE2 residues vs negative RBD residues & vice versa
select ace2_pos, ace2 and resn ARG+LYS+HIS
select ace2_neg, ace2 and resn ASP+GLU

select rbd_pos, rbd and resn ARG+LYS+HIS
select rbd_neg, rbd and resn ASP+GLU

select ace2_salt, byres (ace2_pos within 4 of rbd_neg or ace2_neg within 4 of rbd_pos)
select rbd_salt,  byres (rbd_pos within 4 of ace2_neg or rbd_neg within 4 of ace2_pos)

create 1_4_res, ace2_salt or rbd_salt
dist   1_4_bonds, ace2_salt, rbd_salt, 4.0
show   sticks, 1_4_res
show   dash,   1_4_bonds
color blue, 1_4_res
color blue, 1_4_bonds

# 5. HYDROPHOBIC residues
# Hydrophobic residue list as PyMOL string
set hydro_res, ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO

# 1) Hydrophobic residues of ACE2 contacting RBD
select ace2_hydro, ace2 and resn ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO and byres (ace2 within 4 of rbd)

# 2) Hydrophobic residues of RBD contacting ACE2
select rbd_hydro,  rbd  and resn ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO and byres (rbd  within 4 of ace2)

# 3) Combine all hydrophobic interface residues
create 1_5_hydro_res, ace2_hydro or rbd_hydro

# ---- Visualization ----
show sticks, 1_5_hydro_res
color orange, 1_5_hydro_res
set stick_radius, 0.25, 1_5_hydro_res

# Optional: show hydrophobic–hydrophobic contacts
dist 1_5_hydro_contacts, ace2_hydro, rbd_hydro, 4.0
color orange, 1_5_hydro_contacts
set dash_width, 2, 1_5_hydro_contacts
set dash_gap, 0.25, 1_5_hydro_contacts

# Label hydrophobic residues
label 1_5_hydro_res and name CA, resn + "-" + resi
set label_color, orange, 1_5_hydro_res
set label_size, 14, 1_5_hydro_res

print("Hydrophobic patch residues (ACE2):")
iterate ace2_hydro and name CA, print(resn, resi)

print("Hydrophobic patch residues (RBD):")
iterate rbd_hydro and name CA, print(resn, resi)