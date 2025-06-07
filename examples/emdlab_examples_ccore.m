% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));
% dimensions
ww = 50; % width of the window [mm]
wc = 10; % with of the core [mm]
wb = 10; % width of the bar [mm]
hw = 50; % height of the window [mm]
gap = 0.5; % air gap length [mm]
d0 = 2; % distance between coil arms and C-Core [mm]
wcoil = 10; % width of the coil [mm]
at = 2400; % coil arm ampere turn [A]
ysh = -5; % the amound of the bar shift in y direction [mm]
% geometry generation
g = emdlab_g2d_db;
l1t = g.addClosedPolylineLoop([0,ww+wc,ww+wc,wc,wc,ww+wc,ww+wc,0],[0,0,wc,wc,hw+wc,hw+wc,hw+2*wc,hw+2*wc]);
g.addFace('core', l1t);
l2t = g.addRectangleLoop(wc+d0,wc+d0,wcoil,hw-2*d0);
g.addFace('coil_1', l2t);
l3t = g.addRectangleLoop(-d0-wcoil,wc+d0,wcoil,hw-2*d0);
g.addFace('coil_2', l3t);
l4t = g.addRectangleLoop(wc+ww+gap,ysh,wb,hw+2*wc);
g.addFace('yoke', l4t);
l5t = g.addRectangleLoop(-30,-30,120,130);
g.addFace('air', l5t, l1t, l2t, l3t, l4t);
g.setMeshMaxLength(3);
g.showSketch;
m = g.generateMesh('mg0');
m.setMeshZoneColor('core',90,90,90);
m.setMeshZoneColor('yoke',90,90,90);
m.setMeshZoneColor('air',150,237,239);
m.setMeshZoneColor('coil_1',250,137,39);
m.setMeshZoneColor('coil_2',250,137,39);
% standart refinement of the mesh
m.strefine;
m.strefine;
m.ggmesh;
m.showmzs;
% add material to mesh database
m.addMaterial('air');
m.addMaterial('es_M400_50A');
m.addMaterial('copper');
% set material of the regions
m.setMaterial('air','air');
m.setMaterial('core','es_M400_50A');
m.setMaterial('yoke','es_M400_50A');
m.setMaterial('coil_1','copper');
m.setMaterial('coil_2','copper');
% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
% setting physical units
s.setUnit('length', 'mm');
s.setUnit('currentDensity', 'A/mm^2');
s.setUnit('magneticVectorPotential', 'A/m');
s.setDepth(20);
% define winding
s.defineWinding('winding_1');
s.addMeshZone2Winding('winding_1','coil_1',1,'negative');
s.addMeshZone2Winding('winding_1','coil_2',1,'positive');
s.setWindingCurrent('winding_1', at);
% set boundary condition
s.setAzBC(m.getfbn, 0);
% run solver
figure;
s.setSolverRelativeError(1e-4);
s.solve;
s.plotBmagF(15, 'core', 'yoke');
clim([0,2]);
colormap(jet(15));
