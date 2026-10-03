############################################
# 1. SETUP & CUSTOM COLORS
############################################
reinitialize everything
bg_color white
set ray_trace_mode, 1
set antialias, 2
set ambient, 0.4

# User Defined Palette
set_color l, [0.4980, 0.7882, 0.9373] # Sky Blue (Protein)
set_color y, [0.9765, 0.8667, 0.3020] # Yellow (Crystal Lig)
set_color b, [0.1294, 0.3216, 0.8078] # Dark Blue (Pocket)
set_color o, [0.9373, 0.5176, 0.3255] # Orange (Predict Lig)
set_color p, [0.5700, 0.2500, 0.8000] # Purple (Metal/Ionic)

# Interaction Palette (Harmonious with above)
set_color hb_col,   [0.85, 0.30, 0.30] # Soft Red (H-bonds)
set_color hyd_col,  [0.70, 0.70, 0.70] # Grey (Hydrophobic)
set_color pi_col,   [0.20, 0.60, 0.20] # Forest Green (Pi-Stacking/Cation-Pi)
set_color ion_col,  [0.9765, 0.8667, 0.3020] # Yellow <-- Light Magenta (Ionic)

############################################
# 2. LOAD DATA & DEFINE OBJECTS
############################################
# Load files (Modify filenames if needed)
load ../../Q2_1/1R4L.pdb, 1R4L
load ../../Q2_3/PLIP Input/result.dock4, rs2

# 2.1 Extract Components
# 1R4L (Light Blue)
color l, 1R4L

# Crystal Ligand (Yellow)
create crystal_lig, 1R4L and resn XX5
color y, crystal_lig

# Predicted Ligand (Orange)
create predict_lig, rs, 1, 1
color o, predict_lig

# Zinc (Purple/Blue context)
select zn, 1R4L and resn ZN
color p, zn 

############################################
# 3. DEFINE POCKET (4.5 A)
############################################
# Select protein residues within 4.5A of predicted ligand
select pocket_raw, 1R4L and byres (1R4L within 4.5 of predict_lig and not element H)
create predict_pocket, pocket_raw
# Color Pocket Residues (Dark Blue as requested)
color b, predict_pocket
color b, pocket_raw
show sticks, predict_pocket

############################################
# 4. INTERACTION ANALYSIS (AUTOMATED)
############################################

# 4.1 Hydrogen Bonds (3.4 A, Soft Red)
# mode=2 finds polar contacts
distance dist_hb, predict_lig, predict_pocket, 3.4, mode=2
color hb_col, dist_hb
set dash_width, 3.0, dist_hb
set dash_gap, 0.2, dist_hb

# 4.2 Metal Coordination (Zn, Purple)
distance dist_zn, predict_lig, zn, 4.0
color p, dist_zn
set dash_width, 4.0, dist_zn
set dash_gap, 0.1, dist_zn

# 4.3 Hydrophobic Contacts (4.5 A, Grey)
# Select hydrophobic residues in pocket
select hydro_res, predict_pocket and resn ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO+TYR
distance dist_hyd, predict_lig, hydro_res, 4.5, mode=0
color hyd_col, dist_hyd
set dash_width, 1.0, dist_hyd
set dash_transparency, 0.2, dist_hyd # Semi-transparent background
# Remove H-bond overlaps from hydrophobic view
# (Optional: pure geometric distance usually includes H-bonds, visual separation is key)

# 4.4 Electrostatic Interactions (Salt Bridges)
# Pos/Neg residues within 10.0A
select pos_res, predict_pocket and resn ARG+LYS+HIS
select neg_res, predict_pocket and resn ASP+GLU
distance dist_ion, predict_lig, pos_res or neg_res, 10.0, mode=2
color ion_col, dist_ion
set dash_width, 2.5, dist_ion
set dash_gap, 0.35, dist_ion

# 4.5 Pi-Stacking / Cation-Pi (Aromatic Rings)
# Select aromatic residues (Phe, Tyr, Trp, His) within 5.0A
select arom_res, predict_pocket and resn PHE+TYR+TRP+HIS
# Create a visual distance to potential Pi-interacting residues
distance dist_pi, predict_lig, arom_res, 4.5, mode=0
color pi_col, dist_pi
set dash_width, 2.0, dist_pi
set dash_gap, 0.5, dist_pi

############################################
# 5. VISUALIZATION STYLES
############################################
hide everything, 1R4L and not polymer.protein
show sphere, 1R4L and resn ZN

# Show Protein
show cartoon, protein
set cartoon_transparency, 0.6

# Show Ligands
show sticks, crystal_lig
show sticks, predict_lig
show spheres, zn
set sphere_scale, 0.45, zn

# Show Pocket Residues
show sticks, predict_pocket
set stick_radius, 0.2, predict_pocket
set stick_radius, 0.25, crystal_lig
set stick_radius, 0.25, predict_lig

# Clean up Labels (Hide distance values)
hide labels, dist_*

############################################
# 6. LABELS & LEGEND
############################################
set label_font_id, 7
set label_size, 20
set label_outline_color, white
set label_position, (1, 1, 1)

# Ligand Labels
label crystal_lig and name C1, "Crystal (1R4L)"
label predict_lig and name C1, "MLN-4760-Predicted"

# Residue Labels (Format: HIS-374)
label predict_pocket and name CA, "%s-%s" % (resn, resi)
set label_color, b, predict_pocket

# Legend (Fake atoms to show color coding)
pseudoatom legend_hb, pos=[10,10,0], label="H-Bond"
pseudoatom legend_ele, pos=[10,8,0],  label="Electrostatic Interaction"
pseudoatom legend_pi, pos=[10,6,0],  label="Pi Interaction"
pseudoatom legend_hyd, pos=[10,4,0],  label="Hydrophobic Interaction"
pseudoatom legend_lig, pos=[10,2,0],  label="MLN-4760-Predicted"
pseudoatom legend_lig_exp, pos=[10,2,0],  label="MLN-4760-Experimental"
pseudoatom legend_pro, pos=[10,0,0],  label="ACE2"
pseudoatom legend_zn, pos=[10,12,0],  label="Zn"
set label_color, red, legend_hb
set label_color, y, legend_ele
set label_color, pi_col, legend_pi
set label_color, grey, legend_hyd
set label_color, o, legend_lig
set label_color, y, legend_lig_exp
set label_color, l, legend_pro
set label_color, p, legend_zn


# Final View
zoom predict_lig, 5
clip slab, 10
rebuild