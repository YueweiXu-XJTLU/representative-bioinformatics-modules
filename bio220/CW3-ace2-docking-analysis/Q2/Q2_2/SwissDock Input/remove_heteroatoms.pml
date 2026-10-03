select water, 1R4L_fixed and resn HOH
select NAG_800, 1R4L_fixed and resn NAG and resi 800
select NAG_801, 1R4L_fixed and resn NAG and resi 801
select CL, 1R4L_fixed and resn CL
select ZN, 1R4L_fixed and resn ZN
select XX5, 1R4L_fixed and resn XX5

# protein (no xx5 + zn + water + cl + nag_800 + nag_801)
create 1R4L_protein, 1R4L_fixed and not XX5
save 1R4L_protein.pdb, 1R4L_protein

# protein + zn (no xx5 + water + cl + nag_800 + nag_801)
create 1R4L_protein_zn, 1R4L_fixed and not XX5 and not NAG_800 and not NAG_801 and not CL and not water
save 1R4L_protein_zn.pdb, 1R4L_protein_zn

# protein + zn + water (no xx5 + cl + nag_800 + nag_801)
create 1R4L_protein_zn_water, 1R4L_fixed and not XX5 and not NAG_800 and not NAG_801 and not CL
save 1R4L_protein_zn_water.pdb, 1R4L_protein_zn_water

# protein + zn + water + cl (no xx5 + nag_800 + nag_801)
create 1R4L_protein_zn_water_cl, 1R4L_fixed and not XX5 and not NAG_800 and not NAG_801
save 1R4L_protein_zn_water_cl.pdb, 1R4L_protein_zn_water_cl

# protein + zn + water + cl + nag_800 + nag_801 (no xx5)
create 1R4L_protein_zn_water_cl_nag800_nag801, 1R4L_fixed and not XX5
save 1R4L_protein_zn_water_cl_nag800_nag801.pdb, 1R4L_protein_zn_water_cl_nag800_nag801