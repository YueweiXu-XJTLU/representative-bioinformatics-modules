reinitialize everything

load ../../Q1_1/6M0J.pdb, native_raw
create native, native_raw and polymer.protein
select rec_native, native and chain A and polymer.protein
select lig_native, native and chain E and polymer.protein

load ../cluspro_results/model.000.00.pdb
set_name rec.pdb, rec_dock_raw
set_name lig.000.00.pdb, lig_dock_raw
create dock, rec_dock_raw or lig_dock_raw
select rec_dock, dock and chain A and polymer.protein
select lig_dock, dock and chain E and polymer.protein

# align ACE2 (REC)
align rec_dock and name CA, rec_native and name CA

# Match: read scoring matrix.
# Match: assigning 597 x 597 pairwise scores.
# MatchAlign: aligning residues (597 vs 597)...
# MatchAlign: score 3241.000
# ExecutiveAlign: 597 atoms aligned.
# ExecutiveRMS: 1 atoms rejected during cycle 1 (RMSD=0.00).
# Executive: RMSD =    0.000 (596 to 596 atoms)

# select atoms to align
select native_if_rec, rec_native and resi 24+27+28+30+31+34+35+37+38+41+42+79+82+83+330+353+354+355+357+393
select native_if_lig, lig_native and resi 417+446+449+453+455+456+475+486+487+489+493+496+498+500+501+502+505
create native_if, native_if_rec or native_if_lig
select native_if_CA, native_if and name CA and (alt '' or alt A)

select dock_if_rec, rec_dock and resi 24+27+28+30+31+34+35+37+38+41+42+79+82+83+330+353+354+355+357+393
select dock_if_lig, lig_dock and resi 417+446+449+453+455+456+475+486+487+489+493+496+498+500+501+502+505
create dock_if, dock_if_rec or dock_if_lig
select dock_if_CA, dock_if and name CA

# calculate alpha carbon RMSD
alter native_if_CA, chain='Z'
alter dock_if_CA, chain='Z'
alter all, segi=''
rms_cur (native_if_CA), (dock_if_CA)

# Executive: RMSD =   37.917 (37 to 37 atoms)

# !!!!!
# the rms_* and fit commands required that the atoms in each selection to have matching identifiers: such as chain, resi, segi, resn, etc.
# You can use the "alter" command to set the identifiers to matching values before issuing such commands:
# alter all, segi=''
# alter c321 and chain E, chain='A'
# alter c321 and chain I, chain='B'
# etc.
# !!!!!

############################################
# 1. SHOW & COLOR
############################################
hide everything

show cartoon, native
show cartoon, dock

show sticks, native_if
show sticks, dock_if

# native cartoon (light blue)
set_color a, [0.55, 0.80, 0.95]
# dock cartoon (yellow)
set_color b, [0.98, 0.85, 0.25]
# native interface sticks (deep blue)
set_color c, [0.10, 0.25, 0.80]
# dock interface sticks (orange)
set_color d, [0.95, 0.55, 0.20]

color a, native
color b, dock
color c, native_if
color d, dock_if

############################################
# 2. AESTHETIC SETTINGS
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

set cartoon_transparency, 0.40, native
set cartoon_transparency, 0.40, dock

set label_outline_color, -1
set label_font_id, 7
set label_size,12

label native and name CA and resi 100, "Experimental Complex"
set label_color, a, native

label dock and name CA and resi 451, "ClusPro Model Complex"
set label_color, b, dock

label native_if and name CA and resi 489, "Experimental Interface"
set label_color, c, native_if

label dock_if and name CA and resi 489, "ClusPro Model Interface - Ligand"
set label_color, d, dock_if

label dock_if and name CA and resi 493, "ClusPro Model Interface - Receptor"
set label_color, d, dock_if
