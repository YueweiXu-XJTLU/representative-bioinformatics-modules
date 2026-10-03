reinitialize everything

###########################
# 1. Load structural data
###########################

load 6M0J.pdb, 6m0j
create ace2_rbd, 6m0j and chain A+E and polymer.protein

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

############################################################
# METHOD 2: LITERATURE RESIDUES (Lan et al. 2020)
############################################################

select ace2_lit, ace2 and resi 24+27+28+30+31+34+35+37+38+41+42+79+82+83+330+353+354+355+357+393
select rbd_lit,  rbd  and resi 417+446+449+453+455+456+475+486+487+489+493+496+498+500+501+502+505

color cyan, ace2_lit
color cyan, rbd_lit

# Hydrogen bonds from literature (RBD side)
select rbd_hbond_lit, rbd and resi 417+487+493+505+449+500++501+446+449+489+487+502+505
# Hydrogen-bond partners on ACE2 (from Table 1)
select ace2_hbond_lit, ace2 and resi 24+30+35+37+38+41+42+83+325+329+330+353+393

# Salt bridges (lit: mainly K417)
select rbd_salt_lit, rbd and resi 417
select ace2_salt_lit, ace2 and resi 30+329

# Show sticks for visibility
show sticks, ace2_lit or rbd_lit or rbd_hbond_lit or ace2_hbond_lit or rbd_salt_lit or ace2_salt_lit

# create objects
create 2_res, ace2_lit or rbd_lit
show   sticks, 2_res

create 2_hbond_res, ace2_hbond_lit or rbd_hbond_lit
dist   2_bonds_bonds, ace2_hbond_lit, rbd_hbond_lit, 4.0
show   sticks, 2_hbond_res
show   dash,   2_bonds_bonds

create 2_salt_res, ace2_salt_lit or rbd_salt_lit
dist   2_salt_bonds, ace2_salt_lit, rbd_salt_lit, 4.0
show   sticks, 2_salt_res
show   dash,   2_salt_bonds

# Colouring
color cyan, 2_res
color cyan, 2_hbond_res
color cyan, 2_bonds_bonds
color cyan, 2_salt_res
color cyan, 2_salt_bonds

############################################################
# METHOD 3: PDBePISA-CURATED RESIDUES
############################################################

select ace2_pisa, ace2 and resi 136+150+151+154+155+156+157+158+159+160+246+247+249+250+251+253+279+280+281+283+286+287+295+298+299+601
select rbd_pisa,  rbd  and resi 345+346+349+351+352+357+359+360+394+396+444+445+446+447+448+449+450+452+457+465+466+467+468+469+470+471+472+481+482+483+484+490+492

# PISA hydrogen bonds
select ace2_hbond_pisa, ace2 and resi 150+157+160+281=286+295+299+155+156+247+250
select rbd_hbond_pisa,  rbd  and resi 346+449+444+346+470+468+369+357+450+470+482

# PISA salt bridges
select ace2_salt_pisa, ace2 and resi 150+157+160+295+299
select rbd_salt_pisa,  rbd  and resi 346+357+444

# create objects
create 3_res, ace2_pisa or rbd_pisa
show   sticks, 3_res

create 3_hbond_res, ace2_hbond_pisa or rbd_hbond_pisa
dist   3_bonds_bonds, ace2_hbond_pisa, rbd_hbond_pisa, 4.0
show   sticks, 3_hbond_res
show   dash,   3_bonds_bonds

create 3_salt_res, ace2_salt_pisa or rbd_salt_pisa
dist   3_salt_bonds, ace2_salt_pisa, rbd_salt_pisa, 4.0
show   sticks, 3_salt_res
show   dash,   3_salt_bonds

# Colouring
color violet, 3_res
color violet, 3_hbond_res
color violet, 3_bonds_bonds
color violet, 3_salt_res
color violet, 3_salt_bonds

############################################################
# COMPARATIVE ANALYSIS: CORE + EXTENDED INTERFACE          #
############################################################

# Intersection
select ace2_core, ace2_contact and ace2_bsa and ace2_hbond and ace2_salt and ace2_lit
select rbd_core,  rbd_contact  and rbd_bsa  and rbd_hbond  and rbd_salt  and rbd_lit

color red, ace2_intersection
color red, rbd_intersection

# Union
select ace2_all, ace2_contact or ace2_bsa or ace2_hbond or ace2_salt or ace2_lit
select rbd_all,  rbd_contact  or rbd_bsa  or rbd_hbond  or rbd_salt  or rbd_lit

color green, ace2_union
color green, rbd_union


