%{
note: simulation a cubic magnet in free air
%}

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

% setting the wireframe mesh
g.setLoopMeshMaxLength(l1t,1);
g.setLoopMeshMaxLength(l2t,5);

% mesh generation
m = g.generateMesh('mg0');

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

% apply boundary conditions: zero vector potential on all boundary nodes
s.setAzBC(m.getfbn, 0);

% run solver
s.solve;

% visualize solution
s.gui;