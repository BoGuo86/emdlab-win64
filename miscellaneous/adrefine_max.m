clc
clear

m = TMDBC;
Nm = 8;
glib_spm_rt0(30,40,45,Nm,0.8,5,2);
m.read_g2d_bin('geom.g2d','MG0');

m. ggmesh;
m.showm
