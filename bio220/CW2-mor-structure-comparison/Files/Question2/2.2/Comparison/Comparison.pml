reinitialize everything

############################################
# 1. 加载三个结构
#    - FULL: 8EF6 复合体（实验结构）
#    - AF:   AlphaFold 预测模型（整条，不裁剪）
#    - SM:   sm 模型（整条，不裁剪，编号与8EF6一致）
############################################
load ../../../Question1/8EF6.pdb, FULL
load ../../../Question1/YueweiXu_AlphaFold.pdb, AF
load ../../../Question1/YueweiXu_SwissModel.pdb, SM

############################################
# 2. 从 8EF6 里提取 μOR 主链 (chain R) 作为 RECEPTOR
#    还要暂时提取配体 (MOI) 来界定口袋区域
############################################
create RECEPTOR, FULL and chain R
create LIGAND_TMP, FULL and resn MOI

############################################
# 3. 全局美学参数（白底、柔光）
############################################
bg_color white
set ray_opaque_background, off
set antialias, 4
set ambient, 0.4
set specular, 0.2
set shininess, 10
set reflect, 0
set ray_shadow, off

# 字体样式（跟前面一致）
set label_outline_color, -1
set label_font_id, 7
set label_size, 20

hide everything, all

############################################
# 4. 定义颜色规范
############################################
# 深蓝 = cartoon 主色
set_color col_REC,   [0.1294, 0.3216, 0.8078]
# 浅蓝 = 口袋高亮 / 标签颜色
set_color col_STICK, [0.4980, 0.7882, 0.9373]

############################################
# 5. 把三个结构画成 cartoon，全部着深蓝
#    注意：AF 和 SM 我们不做任何截断，整条都显示
############################################
# RECEPTOR (8EF6链R)
show cartoon, RECEPTOR
set cartoon_smooth_loops, on, RECEPTOR
set cartoon_fancy_helices, on, RECEPTOR
set cartoon_highlight_color, grey90, RECEPTOR
set cartoon_transparency, 0.0, RECEPTOR
color col_REC, RECEPTOR

# AF 全模型
show cartoon, AF
set cartoon_smooth_loops, on, AF
set cartoon_fancy_helices, on, AF
set cartoon_highlight_color, grey90, AF
set cartoon_transparency, 0.0, AF
color col_REC, AF

# SM 全模型
show cartoon, SM
set cartoon_smooth_loops, on, SM
set cartoon_fancy_helices, on, SM
set cartoon_highlight_color, grey90, SM
set cartoon_transparency, 0.0, SM
color col_REC, SM

############################################
# 6. 在 RECEPTOR 上定义口袋区域（几何: 配体4Å内）
#    用于：
#    (a) 让口袋那段 cartoon 半透明 0.5
#    (b) 生成 LABELS_OBJ（口袋关键残基标签）
############################################
# 口袋残基 = 距离配体 ≤4Å 的 RECEPTOR 残基
select POCKET_GEOM_RECEPTOR, byres (RECEPTOR within 4 of LIGAND_TMP)


# 我们构建一个专门的标签对象 LABELS_OBJ，
#   它拷贝口袋残基，以便后面对象列表里能看到它。
create LABELS_OBJ, POCKET_GEOM_RECEPTOR

# 给 LABELS_OBJ 打标签：残基名-编号
# 我们用它的 Cα 来避免每个原子都标
label (LABELS_OBJ and name CA), resn
set label_color, col_STICK, (LABELS_OBJ and name CA)

# 隐藏 LABELS_OBJ 的几何表示，只保留文字标签在图里
hide everything, LABELS_OBJ
# （注意：hide everything 不会自动 hide labels，
#  如果你想以后隐藏标签，用 `hide labels, LABELS_OBJ`）

############################################
# 7. 对齐：用口袋残基编号对齐 AF 和 SM 到 RECEPTOR
#
#    super 的思路：
#    - RECEPTOR / SM / 8EF6 使用同一编号: 149 150 153 238 298 320
#    - AF 使用对应到 8EF6 的映射: 49 50 53 138 198 220
#
#    我们只用这些口袋残基驱动刚体对齐，
#    但 AF / SM 本身保持完整 (整条链不裁剪)。
############################################
# 8EF6 口袋锚定（RECEPTOR上的编号）
select POCKET_E8_ALIGN, RECEPTOR and resi 149+150+153+238+298+320

# AF 口袋锚定（AlphaFold编号，49..220）
select POCKET_AF_ALIGN, AF and resi 49+50+53+138+198+220

# SM 口袋锚定（SM与8EF6一致编号）
select POCKET_SM_ALIGN, SM and resi 149+150+153+238+298+320

# 让 AF 对齐到 RECEPTOR
super AF and POCKET_AF_ALIGN, RECEPTOR and POCKET_E8_ALIGN

# 让 SM 对齐到 RECEPTOR
super SM and POCKET_SM_ALIGN, RECEPTOR and POCKET_E8_ALIGN



show stick, POCKET_AF_ALIGN
show stick, POCKET_SM_ALIGN
set stick_radius, 0.25, POCKET_AF_ALIGN
set stick_radius, 0.25, POCKET_SM_ALIGN
set_color col_STICK, [0.4980, 0.7882, 0.9373]
color col_STICK, POCKET_AF_ALIGN
color col_STICK, POCKET_SM_ALIGN
set cartoon_transparency, 0.5, (AF and POCKET_AF_ALIGN)
set cartoon_transparency, 0.5, (SM and POCKET_SM_ALIGN)
# === 统一并应用标签样式（替换原 label 与 set label_color 段）===
# 基础样式（与全局一致，放在创建标签之前以确保生效）
set label_font_id, 7
set label_size, 20
set label_outline_color, -1

# 创建标签
label (AF and resi 100 and name CA), "AlphaFold"
label (SM and resi 100 and name CA), "SWISS-MODEL"
label (POCKET_AF_ALIGN and name CA), resn
label (POCKET_SM_ALIGN and name CA), resn

# 统一配色：与口袋标签同色（col_STICK）
set label_color, col_REC, (AF and resi 100 and name CA)
set label_color, col_REC, (SM and resi 100 and name CA)
set label_color, col_STICK, POCKET_AF_ALIGN
set label_color, col_STICK, POCKET_SM_ALIGN
set label_color, col_STICK, POCKET_AF_ALIGN
set label_color, col_STICK, POCKET_SM_ALIGN




############################################
# 8. 镜头定位到口袋区域
#    先用配体作为参考点
############################################
center LIGAND_TMP
orient LIGAND_TMP
zoom   LIGAND_TMP, 10

# !!! 统一视角：
# 现在手动旋转到你要的最终出版角度（口袋朝向镜头，三条结构叠好）。
# 然后运行 get_view，复制输出。
# 把 get_view 输出替换掉下面这段占位符。
set_view (\
     0.834626615,   -0.429334015,    0.345066786,\
     0.269538194,    0.864675462,    0.423890591,\
    -0.480361044,   -0.260781854,    0.837402165,\
     0.000000000,    0.000000000, -116.565193176,\
   138.746917725,  125.685546875,  133.688415527,\
  -8042.309570312, 8275.440429688,  -20.000000000 )

############################################
# 11. 出图（开 ray / 导出 png）
############################################
# ray 3000,3000
# png final_AF_SM_RECEPTOR_overlay.png, dpi=300


# surface
create pocket_af, POCKET_AF_ALIGN
set transparency, 0.4, pocket_af
cmd.label("pocket_af","''")
cmd.show_as("surface"   ,"pocket_af")

create pocket_sm, POCKET_SM_ALIGN
set transparency, 0.4, pocket_sm
cmd.label("pocket_sm","''")
cmd.show_as("surface"   ,"pocket_sm")