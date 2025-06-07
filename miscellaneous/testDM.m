clc
clear 
g = GDBC2D;
r = 0.25;
g = g.newcircularcb('c',[0,0],r,5,5);
g = g.newcircularcb('s',[0,0],1,5,5);

g = g.newdDM('s','s','c');

m = MDBCT(g);
m = m.ggmesh;

m.showmesh