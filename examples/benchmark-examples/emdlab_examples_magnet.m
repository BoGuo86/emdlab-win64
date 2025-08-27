% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));
% define geometry data base
g = emdlab_g2d_db;
% add loops
l1t = g.addCenterRectangleLoop(0,0,20,10);
l2t = g.addCenterRectangleLoop(0,0,120,60);
% define faces
g.addFace('Magnet', l1t);
g.addFace('Air', l2t, l1t);
% set maximum mesh lngth
g.setMeshMaxLength(4);
m = g.generateMesh('mg0');
% perform standard mesh refinement
m.strefine;
% set mesh zone colors
m.setMeshZoneColor('Magnet',90,90,90);
m.setMeshZoneColor('Air',150,237,239);
% add materials
m.addMaterial('air', emdlab_mlib_air);
% set mesh zone materials
m.setMaterial('Magnet','air');
m.setMaterial('Air','air');
% define solver
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
% magnetic coercive force, N38UH-80C
Hc = -922100;
% set magnetization
s.setMagnetization('Magnet',Hc,[1,1]);
% set boundary condition
m.ggmesh;
s.setAzBC(m.getfbn, 0);
% runsolver
s.setSolverMaxIteration(100);
s.setSolverRelativeError(1e-5);
s.solve;
s.gui;