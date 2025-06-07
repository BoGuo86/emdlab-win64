% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));
% define geometry data base
g = emdlab_g2d_db;
% add loops
l1 = g.addRectangleLoop(0,0,4,4);
l2 = g.addClosedPolylineLoop([1,3,3,2,2,1], [1,1,3,3,2,2]);
% add faces
g.addFace('Material', l2);
f1 = g.addFace('Region', l1, l2);
f1.addMeshPoints([1,3]);
g.setMeshMaxLength(0.95);
g.showSketch;
m = g.generateMesh('mm');
% perform standard refinment
for i = 1:0, m.strefine; end
% set mesh zonecolor
m.setMeshZoneColor('Material',90,90,90);
m.setMeshZoneColor('Region',150,237,239);
% add materials
m.addMaterial('air');
m.addMaterial('es_M400_50A');
m.addMaterial('iron');
% set material of the regions
m.setMaterial('Region','air');
% m.setMaterial('Material','iron');
m.setMaterial('Material','es_M400_50A');

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
% setting physical units
s.setUnit('length', 'mm');
s.setDepth(1000);
% set boundary conditions
m.showmzs;
m.ggmesh;
s.setAzBC(m.getfbniol_p0u([0,0], [0,1]), 1e-3);
s.setAzBC(m.getfbniol_p0u([4,0], [0,1]), 4e-3);
% run solver
figure;
s.setSolverRelativeError(1e-6);
s.solve;
s.plotBmag;
colormap(jet(10))
clim([0,2]);

s.evalTotalEnergyCoenergy
s.plotAmag;
m.gd2elements;
s2 = emdlab_solvers_ms2d_tl6_ihnl(m);
% setting physical units
s2.setUnit('length', 'mm');
s2.setDepth(1000);
% set boundary conditions
k1 = m.getfbniol_p0u([0,0], [0,1]);
k2 = m.getfbniol_p0u([4,0], [0,1]);
s2.setAzBC(k1, 1e-3);
s2.setAzBC(k2, 4e-3);
% s2.setAzBC(m.getnIndexOnPoint(50,0),0)
s2.setSolverRelativeError(1e-6);
figure
s2.setMonitor(1);
s2.solve;
s2.plotAmag
s2.evalBe

s2.plotBmagSmooth
colormap(jet(10))
clim([0,2]);


