% example information:
% magneto-static analysis of a single-phase transformer

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% variables: dimensions -> [mm], current -> [A]
gv_ww = 50; % width of the window: must be higher than 2*(gv_d0+gv_wc)
gv_hw = 50; % hight of the window: must be higher than 2*gv_d0
gv_wa = 15; % width of the armature
gv_wy = 15; % width of the yoke
gv_wc = 10; % width of the coil
gv_d0 = 1; % distance between coils and core: must be higher than zero
gv_Lstk = 20; % stack length of the core
gv_KENCLOUSER = 3; % this parameter specefies how much enclouser is larger than core in x and y directions 
gv_Np = 200; % number of turns in primary coil
gv_Ip = 1.6; % current of the primary coil
gv_Ns = 150; % number of turns in secondary coil
gv_Is = 0; % current of the secondary coil

% generation of geometry
g = emdlab_g2d_db;

% add loops
l1 = g.addCenterRectangleLoop(0,0,gv_ww,gv_hw);
l2 = g.addCenterRectangleLoop(0,0,gv_ww+2*gv_wa,gv_hw+2*gv_wy);
l3 = g.addCenterRectangleLoop(-gv_ww/2+gv_d0+gv_wc/2,0,gv_wc,gv_hw-2*gv_d0);
l4 = g.addCenterRectangleLoop(-gv_ww/2-gv_d0-gv_wc/2-gv_wa,0,gv_wc,gv_hw-2*gv_d0);
l5 = g.addCenterRectangleLoop(gv_ww/2-gv_d0-gv_wc/2,0,gv_wc,gv_hw-2*gv_d0);
l6 = g.addCenterRectangleLoop(gv_ww/2+gv_d0+gv_wc/2+gv_wa,0,gv_wc,gv_hw-2*gv_d0);
l7 = g.addCenterRectangleLoop(0,0,gv_KENCLOUSER*(gv_ww+2*gv_wa),gv_KENCLOUSER*(gv_hw+2*gv_wy));

% addfaces
g.addFace('PrimaryCoilArm1', l3);
g.addFace('PrimaryCoilArm2', l4);
g.addFace('SecondaryCoilArm1', l5);
g.addFace('SecondaryCoilArm2', l6);
g.addFace('Core', l2, l1);
g.addFace('InAir', l1, l3, l5);
g.addFace('OutAir', l7, l2, l4, l6);

% mesh generation
g.setMeshMaxLength(3);
g.setLoopMeshMaxLength(l7, 8);
m = g.generateMesh('mg0');

% set mesh zone colors
m.setMeshZoneColor('Core',90,90,90);
m.setMeshZoneColor('InAir',0,255,255);
m.setMeshZoneColor('OutAir',0,255,255);
m.setMeshZoneColor('PrimaryCoilArm1',255,137,39);
m.setMeshZoneColor('PrimaryCoilArm2',255,137,39);
m.setMeshZoneColor('SecondaryCoilArm1',255,137,39);
m.setMeshZoneColor('SecondaryCoilArm2',255,137,39);

% add materials
m.addMaterial('m400', emdlab_mlib_es_M400_50A);
m.addMaterial('copper', emdlab_mlib_copper);

% set mesh zone materials
m.setMaterial('Core','m400');
m.setMaterial('PrimaryCoilArm1','copper');
m.setMaterial('PrimaryCoilArm2','copper');
m.setMaterial('SecondaryCoilArm1','copper');
m.setMaterial('SecondaryCoilArm2','copper');

% get the magneto-static solver
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% define coils
s.defineCoil('PrimaryCoil');
s.addMeshZone2Coil('PrimaryCoil','PrimaryCoilArm1',gv_Np,1);
s.addMeshZone2Coil('PrimaryCoil','PrimaryCoilArm2',gv_Np,-1);

s.defineCoil('SecondaryCoil');
s.addMeshZone2Coil('SecondaryCoil','SecondaryCoilArm1',gv_Ns,1);
s.addMeshZone2Coil('SecondaryCoil','SecondaryCoilArm2',gv_Ns,-1);

s.setCoilCurrent('PrimaryCoil', gv_Ip);
s.setCoilCurrent('SecondaryCoil', gv_Is);

% apply boundary conditions

s.setAzBC(m.getfbn, 0);

% solve and plot results
s.setSolverRelativeError(1e-4);
s.setSolverMaxIteration(100);
s.solve;
s.gui;

