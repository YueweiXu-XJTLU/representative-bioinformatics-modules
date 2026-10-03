############################################
# 0. LOAD 1R4L
############################################
reinitialize everything

load 1R4L.pdb, 1R4L
set_name 1R4L, ace2all

############################################
# 1. BASIC SELECTIONS
############################################

# ACE2 protein
select ace2, ace2all and polymer.protein

# Inhibitor MLN-4760 (XX5)
select inhib, ace2all and resn XX5

# Catalytic Zn2+
select zn, ace2all and resn ZN

############################################
# 2. DISTANCE-BASED 4.5 Å POCKET (21 residues)
############################################

select pocket_raw, ace2 and byres (ace2 within 4.5 of inhib)

# Remove Zn and inhibitor from the pocket selection
select pocket, pocket_raw and not (resn XX5 or resn ZN)

############################################
# 3. COLOR DEFINITIONS
############################################

set_color l,   [0.4980, 0.7882, 0.9373]
set_color y,    [0.9765, 0.8667, 0.3020]
set_color b, [0.1294, 0.3216, 0.8078]
set_color o,  [0.9373, 0.5176, 0.3255]
set_color p,    [0.5700, 0.2500, 0.8000]

# Cartoon color
set cartoon_transparency, 0.6
color l, all

# Catalytic Zn coordination residues (dark blue) 
color b, ace2 and resi 374+378+402
color b, zn

# Hydrogen-bond residues (orange)
color o, ace2 and resi 273+345+346+371+375+505+515

# Remaining pocket residues (yellow)
color y, pocket and resi 145+149+344+347+360+361+363+368+504+510+514

# Inhibitor color (light cyan)
color p, inhib

############################################
# 4. DISPLAY STYLES
############################################

hide everything, all

create cartoon_obj, ace2
show cartoon, cartoon_obj

# ACE2 cartoon
show cartoon, ace2

# Show sticks for pocket residues & Zn & inhibitor
show sticks, pocket
show sphere, zn
show sticks, inhib

set stick_radius, 0.20, all

############################################
# 5. HYDROGEN BONDS & VDW CONTACTS
############################################

# 5.1 Hydrogen Bonds (Strict cutoff, Red)
# mode=2 in dist command detects polar contacts only
distance hb_all, pocket, inhib, cutoff=3.2, mode=2
color red, hb_all
set dash_width, 2.5, hb_all
set dash_gap, 0.2, hb_all
hide labels, hb_all

# 5.2 VDW/Hydrophobic Contacts (Loose cutoff, Grey)
# We calculate ALL contacts within 4.0, but style them differently
distance vdw_all, pocket, inhib, cutoff=4.0, mode=0
color grey50, vdw_all
set dash_width, 1.0, vdw_all
set dash_gap, 0.5, vdw_all
hide labels, vdw_all

############################################
# Zn–Inhibitor & Zn–Protein Coordination (full)
############################################
set_color pink,  [0.90, 0.50, 0.90]

# Reduce Zn²⁺ sphere size
set sphere_scale, 0.45, zn
set sphere_color, b, zn

# 1) Zn–Inhibitor (MLN-4760), find closest 2 oxygen atoms
select inhib_ox, inhib and name O*
dist zn_inhib_1, zn, inhib_ox, 4
dist zn_inhib_2, zn, inhib_ox, 4

# 2) Zn–His374 (NE2 usually coordinates)
dist zn_h374, zn, (ace2 and resi 374 and name NE2), 2.8

# 3) Zn–His378 (NE2)
dist zn_h378, zn, (ace2 and resi 378 and name NE2), 2.8

# 4) Zn–Glu402 (OE1 or OE2; auto-detect closest oxygen)
select glu402_o, ace2 and resi 402 and name OE*
dist zn_e402, zn, glu402_o, 2.8

# Styling (deep blue)
color pink, zn_inhib_1
color pink, zn_inhib_2
color pink, zn_h374
color pink, zn_h378
color pink, zn_e402

set dash_width, 3.0, zn_inhib_1
set dash_width, 3.0, zn_inhib_2
set dash_width, 3.0, zn_h374
set dash_width, 3.0, zn_h378
set dash_width, 3.0, zn_e402

hide labels, zn_inhib_1
hide labels, zn_inhib_2
hide labels, zn_h374
hide labels, zn_h378
hide labels, zn_e402

############################################
# 6. LABELS (all color-matched)
############################################

# Global Font Settings
# Font ID 7 is Helvetica Bold (clear for publications)
set label_font_id, 7
set label_size, 20
set label_position, (1, 1, 1)
set label_shadow_mode, 2

# Label ACE2
label ace2 and resi 100 and name CA, "ACE2"
set label_size, 24, ace2 and resi 100 and name CA
set label_color, l, ace2 and resi 100 and name CA

# Label Zn
label zn, "Zn"
set label_color, b, zn

# Label Inhibitor (MLN-4760)
# Only label atom C2 to avoid clutter
label inhib and name C2, "MLN-4760"
set label_color, p, inhib and name C2

# Label Pocket Residues
# Format: RESN-RESI (e.g., HIS-374) on CA atoms
label pocket and name CA, "%s-%s" % (resn, resi)

# Color Labels by Group (using your palette: b, o, y)

# Group 1: Zinc Coordination (Blue)
select label_zn_res, ace2 and resi 374+378+402 and name CA
set label_color, b, label_zn_res

# Group 2: Hydrogen Bond Residues (Orange)
select label_hb_res, ace2 and resi 273+345+346+371+505+515 and name CA
set label_color, o, label_hb_res

# Group 3: Other Pocket Residues (Yellow)
# Select residues in pocket that are NOT in group 1 or 2
select label_other_res, pocket and name CA and not (label_zn_res or label_hb_res)
set label_color, y, label_other_res

# Cleanup temporary selections
delete label_zn_res
delete label_hb_res
delete label_other_res

############################################
# 7. GLOBAL AESTHETICS (High-Res)
############################################
bg_color white
set ray_trace_mode, 1
set antialias, 2

set cartoon_transparency, 0.6