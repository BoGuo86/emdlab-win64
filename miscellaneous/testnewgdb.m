clc
clear
g = GDBC2D;
g = g.newpolygonald('Lshape',inf,[0,0;1,0;2,0;2,1;2,2;1,2;1,1;0,1]);
m = MDBCT(g);clear g;
m = m.ggmesh;
m = m.strefine;
m = m.strefine;
m = m.strefine;
m.showmesh