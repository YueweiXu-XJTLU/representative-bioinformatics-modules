############################################
# 0. REINITIALIZE
############################################
reinitialize everything


############################################
# 1. COLOR DEFINITIONS
############################################

# ACE2 (light blue)
set_color col_ACE2_cartoon, [0.55, 0.80, 0.95]

# RBD (yellow)
set_color col_RBD_cartoon,  [0.98, 0.85, 0.25]

# ACE2 interface sticks (deep blue cartoon)
set_color col_ACE2_interface, [0.10, 0.25, 0.80]

# RBD interface sticks (orange cartoon)
set_color col_RBD_interface,  [0.95, 0.55, 0.20]

# Contact hydrogen bonds (gray)
set_color col_CONTACT_HB, [0.55, 0.55, 0.55]

# Contact salt bridges (black)
set_color col_CONTACT_SB, [0.05, 0.05, 0.05]

# Literature hydrogen bonds (red)
set_color col_LIT_HB, [0.80, 0.15, 0.20]

# Literature salt bridges (purple)
set_color col_LIT_SB, [0.57, 0.25, 0.80]

set label_outline_color, -1
set label_font_id, 7


############################################
# 2. AESTHETIC SETTINGS
############################################
bg_color white
set ray_opaque_background, off
set antialias, 4
set ambient, 0.35
set direct, 0.55
set specular, 0.15
set shininess, 20
set reflect, 0
set ray_shadow, off
set cartoon_transparency, 0.40


############################################
# 3. LOAD STRUCTURE
############################################
load 6M0J.pdb, 6m0j

create ace2_rbd, 6m0j and chain A+E and polymer.protein
create ace2, ace2_rbd and chain A
create rbd,  ace2_rbd and chain E

hide everything
show cartoon, ace2
show cartoon, rbd
color col_ACE2_cartoon, ace2
color col_RBD_cartoon,  rbd

# ACE2 label
label ace2 and name CA and resi 100, "ACE2"
set label_color, col_ACE2_cartoon, ace2
set label_size, 24, ace2

# RBD label
label rbd and name CA and resi 450, "RBD"
set label_color, col_RBD_cartoon, rbd
set label_size, 24, rbd

############################################
# 4. CONTACT INTERFACE (<=4A)
############################################
create ace2_contact, byres (ace2 within 4 of rbd)
create rbd_contact,  byres (rbd  within 4 of ace2)

show sticks, ace2_contact
show sticks, rbd_contact

color col_ACE2_interface, ace2_contact
color col_RBD_interface,  rbd_contact

# Contact interactions
dist contact_hbonds, ace2_contact, rbd_contact, 4
color col_CONTACT_HB, contact_hbonds
set dash_width, 2, contact_hbonds
set dash_gap, 0.25, contact_hbonds
set dash_color, col_CONTACT_HB, contact_hbonds

# Smaller text for distances (contact)
set label_size, 0, contact_hbonds
set label_color, col_CONTACT_HB, contact_hbonds


############################################
# 5. LABEL ALL INTERFACE RESIDUES
############################################
label ace2_contact and name CA, resn + "-" + resi
set label_color, col_ACE2_interface, ace2_contact

label rbd_contact and name CA, resn + "-" + resi
set label_color, col_RBD_interface, rbd_contact


############################################
# 6. LITERATURE HYDROGEN BONDS
############################################
dist hb_01, (chain E and resi 487 and name ND2), (chain A and resi 24 and name OE1)
dist hb_02, (chain E and resi 417 and name NZ),  (chain A and resi 30 and name OD2)
dist hb_03, (chain E and resi 493 and name NE2), (chain A and resi 35 and name OE2)
dist hb_04, (chain E and resi 505 and name OH),  (chain A and resi 37 and name OE2)
dist hb_05, (chain E and resi 449 and name OH),  (chain A and resi 38 and name OD2)
dist hb_06, (chain E and resi 500 and name OG1), (chain A and resi 41 and name OH)
dist hb_07, (chain E and resi 501 and name N),   (chain A and resi 41 and name OH)
dist hb_08, (chain E and resi 446 and name O),   (chain A and resi 42 and name NE2)
dist hb_09, (chain E and resi 449 and name OH),  (chain A and resi 42 and name NE2)
dist hb_10, (chain E and resi 489 and name OH),  (chain A and resi 83 and name OH)
dist hb_11, (chain E and resi 487 and name OD1), (chain A and resi 83 and name OH)
dist hb_12, (chain E and resi 502 and name N),   (chain A and resi 353 and name O)
dist hb_13, (chain E and resi 505 and name OH),  (chain A and resi 393 and name NH2)

group hbonds_lit, hb_01 hb_02 hb_03 hb_04 hb_05 hb_06 hb_07 hb_08 hb_09 hb_10 hb_11 hb_12 hb_13
color col_LIT_HB, hbonds_lit
show dash, hbonds_lit

set dash_color, col_LIT_HB, hbonds_lit
set label_size, 12, hbonds_lit
set label_color, col_LIT_HB, hbonds_lit

# Make lines thicker and clearer
set dash_width, 3, hbonds_lit
set dash_gap, 0.20, hbonds_lit
set dash_color, col_LIT_HB, hbonds_lit

# Larger labels
set label_size, 0, hbonds_lit
set label_color, col_LIT_HB, hbonds_lit

############################################
# 7. LITERATURE SALT BRIDGES
############################################
dist sb_1, (chain E and resi 417 and name NZ), (chain A and resi 30 and name OD1)
dist sb_2, (chain E and resi 417 and name NZ), (chain A and resi 30 and name OD2)

group salt_lit, sb_1 sb_2
color col_LIT_SB, salt_lit
show dash, salt_lit

set dash_color, col_LIT_SB, salt_lit
set label_size, 12, salt_lit
set label_color, col_LIT_SB, salt_lit

set dash_transparency, 0.4, contact_hbonds

# Thicker lines
set dash_width, 3, salt_lit
set dash_gap, 0.20, salt_lit
set dash_color, col_LIT_SB, salt_lit

# Larger labels
set label_size, 0, salt_lit
set label_color, col_LIT_SB, salt_lit

############################################
# 8. FINAL VIEW
############################################
center ace2_rbd
orient ace2_rbd
zoom ace2_rbd, 8

# Export
# ray 4000,4000
# png interface_final.png, dpi=300


# python
# from pymol import cmd, stored
# 
# def print_residues(selection_name):
#     # 1. 初始化存储列表
#     stored.res_list = []
#     
#     # 2. 遍历对象中的所有原子，收集 (残基名, 残基号) 的元组
#     # 使用元组是为了后续方便排序和去重
#     cmd.iterate(selection_name, "stored.res_list.append((resn, resi))")
#     
#     # 3. 去重 (set) 并排序 (sorted)
#     # key=lambda x: int(x[1]) 确保按数字顺序排序 (如 9, 10) 而不是字符串顺序 (10, 9)
#     # 注意：如果残基号包含插入码(如100A)，需移除 int() 转换或使用更复杂的排序逻辑，这里默认是纯数字
#     unique_residues = sorted(list(set(stored.res_list)), key=lambda x: int(x[1]))
#     
#     # 4. 格式化输出 "RESN-RESI"
#     formatted_output = [f"{res[0]}-{res[1]}" for res in unique_residues]
#     
#     # 5. 打印结果
#     print(f"========== {selection_name} ==========")
#     print(f"Total Residues: {len(formatted_output)}")
#     print(", ".join(formatted_output))
#     print("========================================\n")
# 
# # 调用函数分别输出两个对象
# print_residues("ace2_contact")
# print_residues("rbd_contact")
# python end
# ========== ace2_contact ==========
# Total Residues: 20
# GLN-24, THR-27, PHE-28, ASP-30, LYS-31, HIS-34, GLU-35, GLU-37, ASP-38, TYR-41, GLN-42, LEU-79, MET-82, TYR-83, ASN-330, LYS-353, GLY-354, ASP-355, ARG-357, ARG-393
# ========================================
#  
# ========== rbd_contact ==========
# Total Residues: 17
# LYS-417, GLY-446, TYR-449, TYR-453, LEU-455, PHE-456, ALA-475, PHE-486, ASN-487, TYR-489, GLN-493, GLY-496, GLN-498, THR-500, ASN-501, GLY-502, TYR-505
# ========================================