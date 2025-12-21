% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% dimensions & parameters
ww = 50; % width of the window [mm]
wc = 10; % with of the core [mm]
wb = 10; % width of the bar [mm]
hw = 50; % height of the window [mm]
gap = 0.5; % air gap length [mm]
d0 = 2; % distance between coil arms and C-Core [mm]
wcoil = 10; % width of the coil [mm]
at = 2400; % total coil arm ampere turn [A]
ysh = -5; % the amound of the bar shift in y direction [mm]

% geometry generation
g = emdlab_g2d_db;
l1t = g.addClosedPolylineLoop([0,ww+wc,ww+wc,wc,wc,ww+wc,ww+wc,0],[0,0,wc,wc,hw+wc,hw+wc,hw+2*wc,hw+2*wc]);
g.addFace('core', l1t);
l2t = g.addRectangleLoop(wc+d0,wc+d0,wcoil,hw-2*d0);
g.addFace('coilArm_1', l2t);
l3t = g.addRectangleLoop(-d0-wcoil,wc+d0,wcoil,hw-2*d0);
g.addFace('coilArm_2', l3t);
l4t = g.addRectangleLoop(wc+ww+gap,ysh,wb,hw+2*wc);
g.addFace('yoke', l4t);
l5t = g.addRectangleLoop(-30,-30,120,130);
g.addFace('air', l5t, l1t, l2t, l3t, l4t);

% setting the wireframe mesh
g.setMeshMaxLength(3);

% mesh generation
m = g.generateMesh('mg0');
m.setMeshZoneColor('core',90,90,90);
m.setMeshZoneColor('yoke',90,90,90);
m.setMeshZoneColor('air',0,255,255);
m.setMeshZoneColor('coilArm_1',250,137,39);
m.setMeshZoneColor('coilArm_2',250,137,39);

% standart refinement of the mesh
m.strefine;

% add material to mesh database
m.addMaterial('air', emdlab_mlib_air);
m.addMaterial('m400', emdlab_mlib_es_M400_50A);
m.addMaterial('copper', emdlab_mlib_copper);

% set material of the regions
m.setMaterial('air','air');
m.setMaterial('core','m400');
m.setMaterial('yoke','m400');
m.setMaterial('coilArm_1','copper');
m.setMaterial('coilArm_2','copper');

% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(20);

% define winding
s.defineCoil('coil_1');
s.addMeshZone2Coil('coil_1','coilArm_1',1,1);
s.addMeshZone2Coil('coil_1','coilArm_2',1,-1);
s.setCoilCurrent('coil_1', at);

% set boundary condition: zero vector potential on all boundary nodes
s.setAzBC(m.getfbn, 0);

% run solver
s.setSolverRelativeError(1e-4);
s.setSolverMaxIteration(50);
s.solve;

% visualize solution
s.gui;
