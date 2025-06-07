% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));
% define geometry data base
g = emdlab_g2d_db;
% add loops
l1t = g.addRectangleLoop(-10,-5,20,10);
l2t = g.addRectangleLoop(-60,-30,120,60);
% define faces
g.addFace('Magnet', l1t);
g.addFace('Air', l2t, l1t);
% set maximum mesh lngth
g.setMeshMaxLength(4);
g.showSketch;
m = g.generateMesh('mg0');
% perform standard mesh refinement
m.strefine;
m.strefine;
% set mesh zone colors
m.setMeshZoneColor('Magnet',90,90,90);
m.setMeshZoneColor('Air',150,237,239);
m.showmzs;
% add materials
m.addMaterial('air');
% set mesh zone materials
m.setMaterial('Magnet','air');
m.setMaterial('Air','air');
% define solver
s = emdlab_solvers_ms2d_tl3_ihnl(m);
% setting physical units
s.setUnit('length', 'mm');
s.setUnit('currentDensity', 'A/mm^2');
s.setUnit('magneticVectorPotential', 'A/m');
% magnetic coercive force, N38UH-80C
Hc = -922100;
% set magnetization
s.setMagnetization('Magnet',emdlab_solvers_ms2d_magnetization(Hc,[1,1]));
% set boundary condition
m.ggmesh;
s.setAzBC(m.getfbn, 0);
% runsolver
s.setMonitor(0);
s.setSolverMaxIteration(100);
s.setSolverRelativeError(1e-5);
s.solve;
% plot flux lines
s.plotAmag;