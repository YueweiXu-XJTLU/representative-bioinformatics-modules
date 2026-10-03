load YueweiXu_AlphaFold.pdb, AF

hide everything, AF
show cartoon, AF

set cartoon_smooth_loops, on
set cartoon_fancy_helices, on
set cartoon_highlight_color, grey90
set cartoon_transparency, 0.0

bg_color white
set ray_opaque_background, off
set antialias, 4
set ambient, 0.4
set specular, 0.2
set shininess, 10
set reflect, 0
set ray_shadow, off

set_color af_conf_vhigh,      [0.1294, 0.3216, 0.8078]
set_color af_conf_high,       [0.4980, 0.7882, 0.9373]
set_color af_conf_low,        [0.9765, 0.8667, 0.3020]
set_color af_conf_vlow,       [0.9373, 0.5176, 0.3255]

select af_vhigh, AF and b > 90
select af_high,  AF and b < 90.0001 and b > 70
select af_low,   AF and b < 70.0001 and b > 50
select af_vlow,  AF and b < 50.0001

color af_conf_vhigh,     af_vhigh
color af_conf_high,      af_high
color af_conf_low,       af_low
color af_conf_vlow,      af_vlow

orient AF
zoom AF, 30

set label_outline_color, -1
set label_font_id, 7
set label_size, 20
select NTER_CA, (byres (first (AF and name CA))) and name CA
select CTER_CA, (byres (last  (AF and name CA))) and name CA
label NTER_CA, "N"
label CTER_CA, "C"

# ray 3000,3000
# png YueweiXu_AlphaFold.png, dpi=300