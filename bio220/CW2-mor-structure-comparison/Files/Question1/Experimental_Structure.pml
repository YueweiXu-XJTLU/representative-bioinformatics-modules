reinitialize everything

load 8EF6.pdb, FULL

# ------------------------------------------------------------------
# 全局渲染风格（柔光、无硬阴影，跟你SM那张一致）
# ------------------------------------------------------------------
set antialias, 4
set ambient, 0.4
set specular, 0.2
set shininess, 10
set reflect, 0
set ray_shadow, off
bg_color white
set ray_opaque_background, off

# label 外描边去掉（-1表示无描边）
set label_outline_color, -1

hide everything, FULL

# ------------------------------------------------------------------
# 1. 颜色定义
# ------------------------------------------------------------------
set_color col_TM,        [0.1294, 0.3216, 0.8078]   # 深蓝: 7TM bundle
set_color col_PARTNER,   [0.4980, 0.7882, 0.9373]   # 青色: G蛋白+scFv16
set_color col_LOOP,      [0.9765, 0.8667, 0.3020]   # 黄色: loops
set_color col_LIG,       [0.9373, 0.5176, 0.3255]   # 橙色: 小分子

# ------------------------------------------------------------------
# 2. 受体 μOR (链 R)
# ------------------------------------------------------------------
create RECEPTOR, FULL and chain R
show cartoon, RECEPTOR

set cartoon_smooth_loops, on, RECEPTOR
set cartoon_fancy_helices, on, RECEPTOR
set cartoon_highlight_color, grey90, RECEPTOR
set cartoon_transparency, 0.0, RECEPTOR

# 选出跨膜螺旋 vs 其余环区
select TM_HELICES,    RECEPTOR and ss h
select LOOPS_REGIONS, RECEPTOR and not TM_HELICES

# 着色
color col_TM,   TM_HELICES
color col_LOOP, LOOPS_REGIONS

# 注意：我们保持 cartoon，不再用 putty，那种粗细不均的风格已经去掉

# ------------------------------------------------------------------
# 3. G 蛋白 + scFv16
#    Gαi1  (A,F)
#    Gβ1   (B)
#    Gγ2   (C)
#    scFv16 (E)
# ------------------------------------------------------------------
create PARTNERS, FULL and polymer.protein and (chain A+F+B+C+E)
show cartoon, PARTNERS

set cartoon_smooth_loops, on, PARTNERS
set cartoon_fancy_helices, on, PARTNERS
set cartoon_highlight_color, grey90, PARTNERS

color col_PARTNER, PARTNERS
set cartoon_transparency, 0.3, PARTNERS   # 半透明衬背景，不压μOR

# ------------------------------------------------------------------
# 4. 小分子
#    CLR = 胆固醇
#    MOI = orthosteric agonist
# ------------------------------------------------------------------
create LIG_CHOL, FULL and resn CLR
show sticks, LIG_CHOL
color col_LIG, LIG_CHOL
set stick_radius, 0.25, LIG_CHOL
set stick_transparency, 0.3, LIG_CHOL

create LIG_MOI, FULL and resn MOI
show sticks, LIG_MOI
color col_LIG, LIG_MOI
set stick_radius, 0.25, LIG_MOI
set stick_transparency, 0.3, LIG_MOI

# ------------------------------------------------------------------
# 5. 标签（大字号清晰）
#
# 目标：
#  - 字号 = 20（大图里好读）
#  - 颜色跟所属对象一致
#  - 每个 pseudoatom 再沿相机方向 + 平面方向小幅平移，
#    让文本既不被蛋白挡住，又别全挤在一块
# ------------------------------------------------------------------

# μOR
pseudoatom lbl_receptor, RECEPTOR
hide everything, lbl_receptor
show labels, lbl_receptor
label lbl_receptor, "Mu-type opioid receptor (MOR)"
set label_size, 20, lbl_receptor
set label_font_id, 7, lbl_receptor
set label_color, col_TM, lbl_receptor
# 往前(视线方向)推20, 再向左(-15)和上(+10)一点，避免挡配体
translate [-15,10,20], lbl_receptor, camera=1

# Gαi1
pseudoatom lbl_galpha, FULL and chain A+F and polymer.protein
hide everything, lbl_galpha
show labels, lbl_galpha
label lbl_galpha, "Gαi1"
set label_size, 20, lbl_galpha
set label_font_id, 7, lbl_galpha
set label_color, col_PARTNER, lbl_galpha
# 往前推20, 再向右(+15)和上(+5)
translate [15,5,20], lbl_galpha, camera=1

# Gβ1
pseudoatom lbl_gbeta, FULL and chain B and polymer.protein
hide everything, lbl_gbeta
show labels, lbl_gbeta
label lbl_gbeta, "Gβ1"
set label_size, 20, lbl_gbeta
set label_font_id, 7, lbl_gbeta
set label_color, col_PARTNER, lbl_gbeta
# 前+20, 右+20, 下-5
translate [20,-5,20], lbl_gbeta, camera=1

# Gγ2
pseudoatom lbl_ggamma, FULL and chain C and polymer.protein
hide everything, lbl_ggamma
show labels, lbl_ggamma
label lbl_ggamma, "Gγ2"
set label_size, 20, lbl_ggamma
set label_font_id, 7, lbl_ggamma
set label_color, col_PARTNER, lbl_ggamma
# 前+20, 右+20, 上+15 (放在稍高处，跟Gβ1错开)
translate [20,15,20], lbl_ggamma, camera=1

# scFv16
pseudoatom lbl_scfv, FULL and chain E and polymer.protein
hide everything, lbl_scfv
show labels, lbl_scfv
label lbl_scfv, "scFV16"
set label_size, 20, lbl_scfv
set label_font_id, 7, lbl_scfv
set label_color, col_PARTNER, lbl_scfv
# 前+20, 右+5, 上+30 (抬得更高，像头顶挂签)
translate [5,30,20], lbl_scfv, camera=1

# CLR (胆固醇)
pseudoatom lbl_chol, LIG_CHOL
hide everything, lbl_chol
show labels, lbl_chol
label lbl_chol, "CLR"
set label_size, 20, lbl_chol
set label_font_id, 7, lbl_chol
set label_color, col_LIG, lbl_chol
# 前+20, 左-25, 上+25  （把CLR放受体左上角）
translate [-25,25,20], lbl_chol, camera=1

# MOI (激动剂)
pseudoatom lbl_moi, LIG_MOI
hide everything, lbl_moi
show labels, lbl_moi
label lbl_moi, "MOI"
set label_size, 20, lbl_moi
set label_font_id, 7, lbl_moi
set label_color, col_LIG, lbl_moi
# 前+20, 右+10, 下-10 （把MOI放受体右下角），与CLR分开
translate [10,-10,20], lbl_moi, camera=1

# ------------------------------------------------------------------
# 6. 视角 / 拍法
#
# 我给你三个可选预设:
#
#   PRESET_A: 经典 "GPCR横卧在膜里 + Gi在下面" 视角
#   PRESET_B: 斜45°展示细节（适合 figure 里看口袋+Gα插入位点）
#   PRESET_C: 从膜外侧往胞内面看的俯视（展示Gαi1 engagement）
#
# 你只保留一组 turn 即可，其他注释掉。
# ------------------------------------------------------------------

set_view (\
     0.044388786,   -0.952444851,    0.301460862,\
     0.708642483,   -0.182678714,   -0.681509376,\
     0.704169035,    0.243876487,    0.666835248,\
     0.000000000,    0.000000000, -455.299560547,\
   120.052787781,  115.677146912,   95.187400818,\
   358.961639404,  551.637451172,  -20.000000000 )
center RECEPTOR

# ------------------------------------------------------------------
# 7. 高分辨率输出
#    准备交图时再打开下面两行
# ------------------------------------------------------------------
# ray 3000,3000
# png 8EF6.png, dpi=300