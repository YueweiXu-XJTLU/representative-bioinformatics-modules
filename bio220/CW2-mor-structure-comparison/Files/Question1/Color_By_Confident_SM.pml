load YueweiXu_SwissModel.pdb, SM

hide everything, SM
show cartoon, SM

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

set_color sm_conf_vhigh,      [0.1294, 0.3216, 0.8078]
set_color sm_conf_high,       [0.4980, 0.7882, 0.9373]
set_color sm_conf_low,        [0.9765, 0.8667, 0.3020]
set_color sm_conf_vlow,       [0.9373, 0.5176, 0.3255]

select sm_vhigh, SM and b > 0.9
select sm_high,  SM and b < 0.9001 and b > 0.7
select sm_low,   SM and b < 0.7001 and b > 0.5
select sm_vlow,  SM and b < 0.5001

color sm_conf_vhigh,    sm_vhigh
color sm_conf_high,      sm_high
color sm_conf_low,       sm_low
color sm_conf_vlow,      sm_vlow

orient SM
zoom SM, 30
set label_outline_color, -1
set label_font_id, 7
set label_size, 20
select NTER_CA, (byres (first (SM and name CA))) and name CA
select CTER_CA, (byres (last  (SM and name CA))) and name CA
label NTER_CA, "N"
label CTER_CA, "C"

# ray 3000,3000
# png YueweiXu_SwissModel.png, dpi=300