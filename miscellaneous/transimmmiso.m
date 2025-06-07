clc
clear
materialdir = [cd,'\MaterialsData'];

m = TMDBC;
m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'copper');

glib_tl_1ph(100,100,1000, 2000, 40 ,250);
m.read_g2d_bin('geom.g2d', 'MM1');
% m.strefine
% m.strefine
m.ggmesh
% m.showmzs

s = IHNLNRMSTL3(m);
s.setUnit('length', 'mm');
s.setUnit('currentDensity', 'A/mm^2');
s.setUnit('magneticVectorPotential', 'A/m');
% coils of phase a

s.setDepth(1000);

s.defineMatrix('coil', 1, 1);
s.addMeshZone2Matrix('coil', 'c1', 1);
s.addMeshZone2Matrix('coil', 'c2', -1);

s.bcs.clearAllBCs;
s.bcs.setDirichlet(s.m.getfbn,0);

s.setMatrixCurrent('coil', 1);
s.solve

s.evalMatrixLinkageFlux('coil')
close all
m.showmzs
