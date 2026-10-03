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

############################################################
# 1. CONTACT-BASED (<= 4 A)
############################################################
select ace2_contact, byres (ace2 within 4.0 of rbd)
select rbd_contact,  byres (rbd  within 4.0 of ace2)

# interface
create ace2_interface_obj, ace2_contact
create rbd_interface_obj, rbd_contact

create 1_1_res, ace2_contact or rbd_contact
dist   1_1_bonds, ace2_contact, rbd_contact, 4.0
show   sticks, 1_1_res
show   dash,   1_1_bonds
color yellow, 1_1_res
color yellow, 1_1_bonds

############################################################
# 2. BSA-BASED (SASA)  PyMOL approximation of PISA
############################################################
create ace2_monomer, ace2
create rbd_monomer,  rbd

get_area ace2_monomer, load_b=1
get_area rbd_monomer,  load_b=1
get_area ace2, load_b=1
get_area rbd,  load_b=1

python
from pymol import cmd

def get_deltasasa(free_sel, complex_sel):
    free = cmd.get_model(free_sel).atom
    comp = cmd.get_model(complex_sel).atom
    ds = {}
    for a, b in zip(free, comp):
        key = (a.resi, a.resn)
        ds.setdefault(key, 0)
        ds[key] += a.b - b.b  # SASA = SASA_monomer  SASA_complex
    return ds

ds_ace2 = get_deltasasa("ace2_monomer", "ace2")
ds_rbd  = get_deltasasa("rbd_monomer",  "rbd")

ace2_bsa_res = [resi for (resi,resn),v in ds_ace2.items() if v > 1]
rbd_bsa_res  = [resi for (resi,resn),v in ds_rbd.items()  if v > 1]

print("ACE2 SASA residues:", ace2_bsa_res)
print("RBD  SASA residues:", rbd_bsa_res)

cmd.select("ace2_bsa", "ace2 and resi " + "+".join(ace2_bsa_res))
cmd.select("rbd_bsa",  "rbd and resi "  + "+".join(rbd_bsa_res))
python end

create 1_2_res, ace2_bsa or rbd_bsa
show   sticks, 1_2_res
color green, 1_2_res

############################################################
# 3. HYDROGEN BONDS
############################################################
dist hbonds, ace2, rbd, mode=2

select ace2_hbond, byres (ace2 within 3.2 of rbd)
select rbd_hbond,  byres (rbd  within 3.2 of ace2)

create 1_3_res, ace2_hbond or rbd_hbond
dist   1_3_bonds, ace2_hbond, rbd_hbond, 3.2
show   sticks, 1_3_res
show   dash,   1_3_bonds
color red, 1_3_res
color red, 1_3_bonds

############################################################
# 4. ELECTROSTATIC SALT BRIDGES (<= 10 A)
############################################################
select ace2_pos, ace2 and resn ARG+LYS+HIS
select ace2_neg, ace2 and resn ASP+GLU

select rbd_pos, rbd and resn ARG+LYS+HIS
select rbd_neg, rbd and resn ASP+GLU

select ace2_salt, byres (ace2_pos within 10 of rbd_neg or ace2_neg within 4 of rbd_pos)
select rbd_salt,  byres (rbd_pos within 10 of ace2_neg or rbd_neg within 4 of ace2_pos)

create 1_4_res, ace2_salt or rbd_salt
dist   1_4_bonds, ace2_salt, rbd_salt, 10.0
show   sticks, 1_4_res
show   dash,   1_4_bonds
color blue, 1_4_res
color blue, 1_4_bonds

############################################
# 6. AESTHETIC ADJUSTMENTS (ONLY VISUAL)
############################################

bg_color white
set ray_opaque_background, off
set antialias, 4
set ambient, 0.35
set direct, 0.55
set specular, 0.20
set shininess, 20
set reflect, 0
set ray_shadow, off

set cartoon_transparency, 0.40, ace2_rbd
set transparency, 0.60, 1_2_res

set transparency, 0.4, 1_2_res

set label_outline_color, -1
set label_font_id, 7

show sticks, ace2_interface_obj
show sticks, rbd_interface_obj
show sticks, 1_1_res
show surface, 1_2_res
show sticks, 1_3_res
show sticks, 1_4_res
show sticks, 1_5_hydro_res


set_color col_ACE2_cartoon,   [0.4980, 0.7882, 0.9373]
set_color col_RBD_cartoon,    [0.9765, 0.8667, 0.3020]
set_color col_ACE2_interface, [0.1294, 0.3216, 0.8078]
set_color col_RBD_interface,  [0.9373, 0.5176, 0.3255]

set_color col_purple,    [0.5700, 0.2500, 0.8000]
set_color col_green,     [0.2000, 0.7000, 0.3500]
set_color col_raspberry, [0.9000, 0.3000, 0.5000]
set_color col_extra,     [0.0000, 0.6500, 0.7000]

set_color col_CONTACT_BOND, [0.55, 0.55, 0.55]


color col_ACE2_cartoon, ace2
color col_RBD_cartoon,  rbd

color col_ACE2_interface, ace2_interface_obj
color col_RBD_interface,  rbd_interface_obj

color col_CONTACT_BOND, 1_1_bonds
set dash_color, col_CONTACT_BOND, 1_1_bonds

color col_ACE2_interface, 1_2_res and chain A
color col_RBD_interface, 1_2_res and chain E

color red, 1_3_bonds
set dash_color, red, 1_3_bonds

color blue, 1_4_bonds
set dash_color, blue, 1_4_bonds

set dash_width, 3, 1_3_bonds
set dash_width, 3, 1_5_hydro_contacts

############################################
# LEGEND LABELS (CLEAN VERSION)
############################################

# 1) ACE2
create lab_ace2, ace2 and name CA and resi 100
label lab_ace2, "ACE2"
set label_color, col_ACE2_cartoon, lab_ace2
set label_size, 24, lab_ace2

# 2) RBD
create lab_rbd, rbd and name CA and resi 450
label lab_rbd, "RBD"
set label_color, col_RBD_cartoon, lab_rbd
set label_size, 24, lab_rbd

# 3) HBonds
set label_size, 0, 1_1_bonds
create lab_hbond, first (1_3_res and name CA)
label lab_hbond, "H-Bonds"
set label_color, red, lab_hbond
set label_size, 18, lab_hbond
set label_size, 0, 1_1_bonds

# 4) VDW Interaction  (no label on gray dashed lines)
set label_size, 0, 1_3_bonds
create lab_vdw, first (1_1_res and name CA)
label lab_vdw, "VDW Interaction"
set label_color, grey, lab_vdw
set label_size, 18, lab_vdw

# 5) BSA Contributed Atoms - ACE2
create lab_bsa_ace2, first (1_2_res and chain A and name CA)
label lab_bsa_ace2, "BSA Contributed Atoms - ACE2"
set label_color, col_ACE2_interface, lab_bsa_ace2
set label_size, 18, lab_bsa_ace2

# 6) BSA Contributed Atoms - RBD
create lab_bsa_rbd, first (1_2_res and chain E and name CA)
label lab_bsa_rbd, "BSA Contributed Atoms - RBD"
set label_color, col_RBD_interface, lab_bsa_rbd
set label_size, 18, lab_bsa_rbd

############################################
# 6.6 FINAL VIEW
############################################
center ace2_rbd
orient ace2_rbd
zoom ace2_rbd, 8

# ray 4000,4000
# png cluspro_interface_all_interactions.png, dpi=300