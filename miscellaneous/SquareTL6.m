clc
clear
materialdir = [cd,'\MaterialsData'];
g = GDBC2D;

g = g.newpolygonald('s',inf,[0,0;1,0;1,1;0,1]);

m = MDBCT(g);clear g;
m = m.addMaterial(materialdir,'air');
s = IHLTHSTL6(m);
s.m = s.m.ggmesh;
s.m = s.m.strefine;
s.m = s.m.strefine;
s.m = s.m.strefine;
s.m = s.m.strefine;

s.m = s.m.gd2elements;

k1 = s.m.getIndexOnRay([0,0],[1,0]);
k2 = s.m.getIndexOnRay([1,0],[1,1]);
k3 = s.m.getIndexOnRay([1,1],[0,1]);
k4 = s.m.getIndexOnRay([0,0],[0,1]);

bc = @(p) sin(pi*p(:,2));
    s = s.setdbc(k2,0);
    bc = @(p) sin(pi*p(:,1)); 
    s = s.setdbc(k3,0);
    bc = @(p) sin(pi*p(:,2));
    s = s.setdbc(k4,0);
    bc = @(p) sin(pi*p(:,1)); 
    s = s.setdbc(k1,0);
    s = s.setExcitation('s',10);
s.m = s.m.evalKeFeC('TL6');

s = s.assignEdata;
    s = s.solve;
    s.plotQmag
    
    s = s.evalEtot;
    s.Etot