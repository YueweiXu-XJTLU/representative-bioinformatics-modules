reinitialize everything

# ---------- 1. 加载并重命名为干净对象名 ----------
load 8EF6.A.cif, EXP_8EF6
load AF_AFP35372F1.A.cif, AF_AlphaFold_aligned_to_8EF6
load YueweiXu_SwissModel.pdb.A.cif, SM_SWISS-MODEL_aligned_to_8EF6

set_name EXP_8EF6, EXP
set_name AF_AlphaFold_aligned_to_8EF6, AF
set_name SM_SWISS-MODEL_aligned_to_8EF6, SM

# ---------- 2. 全局渲染风格 ----------
bg_color white
set ray_opaque_background, off
set antialias, 4
set ambient, 0.4
set specular, 0.2
set shininess, 10
set reflect, 0
set ray_shadow, off

# ---------- 3. 只显示 cartoon ----------
hide everything, all

show cartoon, EXP
show cartoon, AF
show cartoon, SM

set cartoon_smooth_loops, on
set cartoon_fancy_helices, on
set cartoon_highlight_color, grey90

# ---------- 4. 定义颜色并上色 ----------
# 深蓝  -> EXP (实验结构)
# 青蓝  -> AF  (AlphaFold)
# 黄色  -> SM  (SWISS-MODEL)
set_color col_EXP, [0.1294, 0.3216, 0.8078]
set_color col_AF,  [0.4980, 0.7882, 0.9373]
set_color col_SM,  [0.9765, 0.8667, 0.3020]

color col_EXP, EXP
color col_AF,  AF
color col_SM,  SM

# ---------- 5. 对齐区段 vs 非对齐区段，分别调透明度 ----------
# 你的对齐信息（8EF6编号）:
# AF aligned: 65-263,267-269,270-349,353
# SM aligned: 65-349,353
# 对 EXP，我们用 SM 的主对齐范围 65-349,353 作为核心显示

# EXP
select EXP_aligned,   EXP and chain R and (resi 65-349+353)
select EXP_unaligned, EXP and not EXP_aligned
set cartoon_transparency, 0.0, EXP_aligned
set cartoon_transparency, 0.5, EXP_unaligned

# AF
select AF_aligned,   AF and (resi 65-263+267-269+270-349+353)
select AF_unaligned, AF and not AF_aligned
set cartoon_transparency, 0.0, AF_aligned
set cartoon_transparency, 0.5, AF_unaligned

# SM
select SM_aligned,   SM and (resi 65-349+353)
select SM_unaligned, SM and not SM_aligned
set cartoon_transparency, 0.0, SM_aligned
set cartoon_transparency, 0.5, SM_unaligned

# ---------- 6. 把相机围绕实验结构来摆好基本视角 ----------
center EXP
orient EXP
# 注意：不要 zoom 这里；我们会在 ray 前手动 zoom

# ---------- 7. Label，和你8EF6那张保持一致的风格 ----------
set label_outline_color, -1
set label_font_id, 7
set label_size, 20

# 为了在屏幕上放三行标签，我们用伪原子作为锚点，然后把锚点往屏幕四角平移
# 颜色 = 该模型主色

pseudoatom lbl_exp
show labels, lbl_exp
label lbl_exp, "Experimental 8EF6"
set label_color, col_EXP, lbl_exp
translate [50,-25,-50], lbl_exp, camera=1

pseudoatom lbl_af
show labels, lbl_af
label lbl_af, "AlphaFold"
set label_color, col_AF, lbl_af
translate [50,-45,-50], lbl_af, camera=1

pseudoatom lbl_sm
show labels, lbl_sm
label lbl_sm, "SWISS-MODEL"
set label_color, col_SM, lbl_sm
translate [50,-35, -50], lbl_sm, camera=1

# ---------- 8. 出图前最后两步（手动执行） ----------
set_view (\
     0.044388786,   -0.952444851,    0.301460862,\
     0.708642483,   -0.182678714,   -0.681509376,\
     0.704169035,    0.243876487,    0.666835248,\
     0.000000000,    0.000000000, -455.299560547,\
   137.479080200,  120.369415283,  125.026123047,\
   358.961669922,  551.637451172,  -20.000000000 )
center EXP and chain R
# 当你把角度手动转到满意的构图后，再执行下面两行：
ray 3000,3000
png superposition.png, dpi=300
save YueweiXu_Alignment_Session.pse