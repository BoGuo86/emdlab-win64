%% initialization
clear all
clc
materialdir = [cd,'\materials'];
%% Mesh
m = TMDBC;
% adding model materials
m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'m19_24ga');
m.addMaterial(materialdir,'iron');
m.addMaterial(materialdir,'laminatedIron');
m.addMaterial(materialdir,'copper');
% reading mesh file
m.read_msh_bin('C:\Users\AliJamalifard\Desktop\gmesh\c_core_s.msh');
% specifing mesh zone materials
m.setMaterial('Zone1','m19_24ga');
% m.setMaterial('Zone1','laminatedIron');
m.setMaterial('Zone2','copper');
m.setMaterial('Zone4','copper');
%% Getting Solver
s = IHNLNRMDWTMTL3(m);
% s = IHLMDWTMTL3(m);
s.setUnit('length', 'mm');
s.setUnit('currentDensity', 'A/mm^2');
s.setUnit('magneticVectorPotential', 'A/m');
s.setUnit('time', 'ms');
s.setDepth(1000);
%% Boundary condition
s.m.strefine;
s.m.ggmesh;
s.bcs.clearAllBCs;
s.bcs.setDirichlet(s.m.getfbn,0);
%% Setting Excitations and solve
s.defineCoil('Zone2', 'turns', 1);
s.defineCoil('Zone4', 'turns', 1, 'direction', 'negative');
ppp = 1;
tmp = linspace(0,ppp,500);
tmp = [tmp;10*sign(sin(2*pi*(1/ppp)*tmp))];
vc = tf_xy('xyData', tmp');
s.defineWinding('win', 'exValue', @(t) sin(2*pi*50*t), 'exType', 'voltage', 'wType', 'solid');
s.addCoil2Winding('win', 'Zone2');
s.addCoil2Winding('win', 'Zone4');
s.setSimulationTime(1, 200);

s.saveBe(1)
s.solve

