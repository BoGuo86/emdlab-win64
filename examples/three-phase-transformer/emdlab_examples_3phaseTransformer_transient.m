% example information:
% transient analysis of a three-phase transformer
% primary coils are connected to the current source
% secondary coils are short circuit -> short circuit test

% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

% variables: dimensions -> [mm], current -> [A]
gv_ww = 835; % width of the window: must be higher than 2*(gv_d0+gv_wc)
gv_hw = 1635; % hight of the window: must be higher than 2*gv_d0
gv_wa = 664; % width of the armature
gv_wy = 664; % width of the yoke
gv_wc = 100; % width of the coil: is the same for primary and secondary
gv_d0 = 67; % distance between coils and core: must be higher than zero
gv_d1 = 50; % distance between coils and core: must be higher than zero
gv_d2 = 267; % distance between coils and core: must be higher than zero
gv_Lstk = 664; % stack length of the core
gv_KENCLOUSER = 2; % this parameter specefies how much enclouser is larger than core in x and y directions 
gv_Np = 200; % number of turns in primary coil
gv_Vp = 2000; % peak current of primary coils
gv_freq = 50; % excitation frequency
gv_Ns = 20; % number of turns in secondary coil

% generation of geometry
g = emdlab_g2d_db;

% add loops
l1 = g.addCenterRectangleLoop(0,0,2*gv_ww+3*gv_wa,gv_hw+2*gv_wy);
l2 = g.addCenterRectangleLoop(-gv_ww/2-gv_wa/2,0,gv_ww,gv_hw);
l3 = g.addCenterRectangleLoop(gv_ww/2+gv_wa/2,0,gv_ww,gv_hw);

% draw coils loops: we have 12 loops
x_centers = -gv_wa/2-gv_d0-gv_d1-1.5*gv_wc;
x_centers(2) = x_centers(end) + gv_d1 + gv_wc;
x_centers(3:4) = -[x_centers(2),x_centers(1)];

lt_coils = zeros(1,12);
% first leg
x_tmp = x_centers - gv_ww - gv_wa;
for i = 1:4
    lt_coils(i) = g.addCenterRectangleLoop(x_tmp(i),0,gv_wc,gv_hw-2*gv_d2);
end
% second leg
x_tmp = x_centers;
for i = 1:4
    lt_coils(i+4) = g.addCenterRectangleLoop(x_tmp(i),0,gv_wc,gv_hw-2*gv_d2);
end
% third leg
x_tmp = x_centers + gv_ww + gv_wa;
for i = 1:4
    lt_coils(i+8) = g.addCenterRectangleLoop(x_tmp(i),0,gv_wc,gv_hw-2*gv_d2);
end

l4 = g.addCenterRectangleLoop(0,0,gv_KENCLOUSER*(2*gv_ww+3*gv_wa),gv_KENCLOUSER*(gv_hw+2*gv_wy));

% add faces
% first leg
g.addFace('PrimaryCoilArm1_A', lt_coils(1));
g.addFace('SecondaryCoilArm1_A', lt_coils(2));
g.addFace('SecondaryCoilArm2_A', lt_coils(3));
g.addFace('PrimaryCoilArm2_A', lt_coils(4));

% second leg
g.addFace('PrimaryCoilArm1_B', lt_coils(1+4));
g.addFace('SecondaryCoilArm1_B', lt_coils(2+4));
g.addFace('SecondaryCoilArm2_B', lt_coils(3+4));
g.addFace('PrimaryCoilArm2_B', lt_coils(4+4));

% third leg
g.addFace('PrimaryCoilArm1_C', lt_coils(1+8));
g.addFace('SecondaryCoilArm1_C', lt_coils(2+8));
g.addFace('SecondaryCoilArm2_C', lt_coils(3+8));
g.addFace('PrimaryCoilArm2_C', lt_coils(4+8));

g.addFace('Core', l1, l2, l3);
g.addFace('InAir1', l2, lt_coils(3:6));
g.addFace('InAir2', l3, lt_coils(7:10));
g.addFace('OutAir', l4, l1, lt_coils([1,2,11,12]));

% mesh generation
g.setMeshMaxLength(gv_wa/10);
g.setLoopMeshMaxLength(l4, 300);
m = g.generateMesh('mg0');

% set mesh zone colors
m.setMeshZoneColor('Core',90,90,90);
m.setMeshZoneColor('InAir1',0,255,255);
m.setMeshZoneColor('InAir2',0,255,255);
m.setMeshZoneColor('OutAir',0,255,255);
m.setMeshZoneColor('PrimaryCoilArm1_A',255,137,39);
m.setMeshZoneColor('PrimaryCoilArm2_A',255,137,39);
m.setMeshZoneColor('SecondaryCoilArm1_A',242,87,95);
m.setMeshZoneColor('SecondaryCoilArm2_A',242,87,95);
m.setMeshZoneColor('PrimaryCoilArm1_B',255,137,39);
m.setMeshZoneColor('PrimaryCoilArm2_B',255,137,39);
m.setMeshZoneColor('SecondaryCoilArm1_B',242,87,95);
m.setMeshZoneColor('SecondaryCoilArm2_B',242,87,95);
m.setMeshZoneColor('PrimaryCoilArm1_C',255,137,39);
m.setMeshZoneColor('PrimaryCoilArm2_C',255,137,39);
m.setMeshZoneColor('SecondaryCoilArm1_C',242,87,95);
m.setMeshZoneColor('SecondaryCoilArm2_C',242,87,95);

% add materials
m.addMaterial('iron', emdlab_mlib_iron);
m.addMaterial('copper', emdlab_mlib_copper);

% set mesh zone materials
m.setMaterial('Core','iron');
m.setMaterial('PrimaryCoilArm1_A','copper');
m.setMaterial('PrimaryCoilArm2_A','copper');
m.setMaterial('SecondaryCoilArm1_A','copper');
m.setMaterial('SecondaryCoilArm2_A','copper');
m.setMaterial('PrimaryCoilArm1_B','copper');
m.setMaterial('PrimaryCoilArm2_B','copper');
m.setMaterial('SecondaryCoilArm1_B','copper');
m.setMaterial('SecondaryCoilArm2_B','copper');
m.setMaterial('PrimaryCoilArm1_C','copper');
m.setMaterial('PrimaryCoilArm2_C','copper');
m.setMaterial('SecondaryCoilArm1_C','copper');
m.setMaterial('SecondaryCoilArm2_C','copper');

% get the magneto-transient solver -> without motion
s = emdlab_solvers_mt2d_tl3_ihnlwtm(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);

% define coils -> fill factor is 0.5
s.defineCoil('PrimaryCoil_A', 'voltage', 'stranded');
s.addMeshZone2Coil('PrimaryCoil_A','PrimaryCoilArm1_A',gv_Np,1,0.5);
s.addMeshZone2Coil('PrimaryCoil_A','PrimaryCoilArm2_A',gv_Np,-1,0.5);
s.defineCoil('SecondaryCoil_A', 'voltage', 'stranded');
s.addMeshZone2Coil('SecondaryCoil_A','SecondaryCoilArm1_A',gv_Ns,-1,0.5);
s.addMeshZone2Coil('SecondaryCoil_A','SecondaryCoilArm2_A',gv_Ns,1,0.5);

s.defineCoil('PrimaryCoil_B', 'voltage', 'stranded');
s.addMeshZone2Coil('PrimaryCoil_B','PrimaryCoilArm1_B',gv_Np,1,0.5);
s.addMeshZone2Coil('PrimaryCoil_B','PrimaryCoilArm2_B',gv_Np,-1,0.5);
s.defineCoil('SecondaryCoil_B', 'voltage', 'stranded');
s.addMeshZone2Coil('SecondaryCoil_B','SecondaryCoilArm1_B',gv_Ns,-1,0.5);
s.addMeshZone2Coil('SecondaryCoil_B','SecondaryCoilArm2_B',gv_Ns,1,0.5);

s.defineCoil('PrimaryCoil_C', 'voltage', 'stranded');
s.addMeshZone2Coil('PrimaryCoil_C','PrimaryCoilArm1_C',gv_Np,1,0.5);
s.addMeshZone2Coil('PrimaryCoil_C','PrimaryCoilArm2_C',gv_Np,-1,0.5);
s.defineCoil('SecondaryCoil_C', 'voltage', 'stranded');
s.addMeshZone2Coil('SecondaryCoil_C','SecondaryCoilArm1_C',gv_Ns,-1,0.5);
s.addMeshZone2Coil('SecondaryCoil_C','SecondaryCoilArm2_C',gv_Ns,1,0.5);

s.setCoilVoltage('PrimaryCoil_A',@(t) gv_Vp*cos(2*pi*gv_freq*t));
s.setCoilVoltage('PrimaryCoil_B',@(t) gv_Vp*cos(2*pi*gv_freq*t-2*pi/3));
s.setCoilVoltage('PrimaryCoil_C',@(t) gv_Vp*cos(2*pi*gv_freq*t+2*pi/3));

s.defineStarConnection('y1', 'PrimaryCoil_' + ["A", "B", "C"])

% apply boundary conditions
s.setAzBC(m.getfbn, 0);

% solve for 800ms, time step is 1ms
s.setMonitor(0);
s.solve(800e-3,1e-3);

g.showSketch;
m.showmzs;
s.plotCoilCurrents
s.plotCoilFluxLinkages

