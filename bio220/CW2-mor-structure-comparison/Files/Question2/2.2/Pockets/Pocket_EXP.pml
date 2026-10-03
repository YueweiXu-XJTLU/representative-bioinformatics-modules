reinitialize everything

############################################
# 1. 加载原始结构
############################################
load ../../../Question1/8EF6.pdb, FULL

############################################
# 2. 单独定义受体(μOR)和配体(MOI)
#    受体：只要蛋白主链，不要小分子
#    μOR 在 8EF6 里通常是链 R 或 M；先试 R
############################################
create RECEPTOR, FULL and chain R and polymer.protein
create BACKGROUND, FULL and not resn MOI and not resn CLR and not chain R

# 配体（正构激动剂）
create LIGAND, FULL and resn MOI

############################################
# 3. 全局绘图风格（柔和白底）
############################################
bg_color white
set ray_opaque_background, off
set antialias, 4
set ambient, 0.4
set specular, 0.2
set shininess, 10
set reflect, 0
set ray_shadow, off

set label_outline_color, -1
set label_font_id, 7
set label_size, 20

hide everything, all

############################################
# 4. 画μOR主结构成 cartoon（深蓝）
############################################
show cartoon, RECEPTOR
set cartoon_smooth_loops, on, RECEPTOR
set cartoon_fancy_helices, on, RECEPTOR
set cartoon_highlight_color, grey90, RECEPTOR
set cartoon_transparency, 0.0, RECEPTOR

set_color col_REC, [0.1294, 0.3216, 0.8078]   # 深蓝
color col_REC, RECEPTOR

show cartoon, BACKGROUND
set cartoon_smooth_loops, on, BACKGROUND
set cartoon_fancy_helices, on, BACKGROUND
set cartoon_highlight_color, grey90, BACKGROUND
set cartoon_transparency, 0.60, BACKGROUND
color col_REC, BACKGROUND

############################################
# 5. 画配体（MOI）成 sticks（橙色）
############################################
show sticks, LIGAND
set stick_radius, 0.25, LIGAND
set stick_transparency, 0.0, LIGAND

set_color col_LIG, [0.9373, 0.5176, 0.3255]   # 橙色
color col_LIG, LIGAND

############################################
# 6. 定义结合位点残基 SITE
#    思路：找出所有在“配体4 Å范围内”的全残基
#    然后只保留蛋白部分（避免把配体自己、别的脂质一起拉进来）
############################################
# 第一步：从配体出发，抓住所有距离 ≤4 Å 的原子，并扩展到整条残基
create SITE, byres LIGAND around 4 and RECEPTOR and chain R

# 对 RECEPTOR 内这些口袋残基对应的 cartoon 提高透明度
select POCKET_RES, byres LIGAND around 4 and RECEPTOR and chain R
set cartoon_transparency, 0.5, (RECEPTOR and POCKET_RES)

############################################
# 7. 把 SITE 这些关键残基高亮成 sticks（也是深蓝）
#    这些就是“正构口袋里直接接触配体的关键氨基酸”
############################################
show sticks, SITE
hide cartoon, SITE   
set stick_radius, 0.25, SITE
set_color col_STICK, [0.4980, 0.7882, 0.9373]
color col_STICK, SITE

############################################
# 8. 加标签（配体 + 口袋残基）
#    - 给配体打 "MOI"
#    - 给口袋残基打残基名+编号（然后你可以手动删掉多余的，只保留2~3个最关键的）
############################################

set label_outline_color, -1
set label_font_id, 7
set label_size, 20
set label_color, col_LIG

# 配体标签
label LIGAND and name C1, "MOI"

# 口袋残基标签
label SITE and name CA, resn + "-" + resi
set label_color, col_STICK, SITE

############################################
# 8.5 画 pocket 内 LIGAND 和 SITE 的近程相互作用
#     目的：显示潜在结合键/接触位点
############################################

# 可选但推荐：先加氢，提高氢键识别准确度
# （如果加了氢导致画面太乱，你也可以注释掉这一行）
# h_add LIGAND
# h_add SITE

# 生成距离/接触对象：CONTACTS
# 解释：
#   distance CONTACTS, <选择1>, <选择2>, 3.6, mode=2
#   - cutoff 3.6 Å （你可以换成 4.0，如果想更松一点）
#   - mode=2 表示只画 LIGAND 和 SITE 之间的相互作用，而不画各自内部的
distance CONTACTS, LIGAND, SITE, 4.0, mode=2

# 美化这些虚线，让它们看起来像“可能的氢键/疏水接触”
set dash_width, 2, CONTACTS
set dash_gap, 0.2, CONTACTS
set dash_color, grey50, CONTACTS

# 调整距离数字标签的可读性
set label_font_id, 7, CONTACTS
set label_size, 14, CONTACTS
set label_outline_color, -1, CONTACTS
set label_color, grey50, CONTACTS


############################################
# 9. 聚焦、放大口袋
############################################
center LIGAND
orient LIGAND
zoom LIGAND, 10


set_view (\
     0.834626615,   -0.429334015,    0.345066786,\
     0.269538194,    0.864675462,    0.423890591,\
    -0.480361044,   -0.260781854,    0.837402165,\
     0.000000000,    0.000000000, -116.565193176,\
   138.746917725,  125.685546875,  133.688415527,\
  -8042.309570312, 8275.440429688,  -20.000000000 )


############################################
# 10. 最终出图（准备导出时再开）
############################################
# ray 3000,3000
# png binding_pocket_exp.png, dpi=300


# surface
create pocket_e8, SITE
set transparency,0.6, pocket_E8
cmd.label("pocket_e8","''")
cmd.show_as("surface"   ,"pocket_e8")