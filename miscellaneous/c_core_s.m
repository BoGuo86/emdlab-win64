%% initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
%% Mesh
m = TMDBC;
m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'laminatedIron');
m.addMaterial(materialdir,'copper');

m.read_msh_bin('C:\Users\AliJamalifard\Desktop\gmesh\c_core_s.msh');

m.setMaterial('Zone1','laminatedIron');
m.setMaterial('Zone2','copper');
m.setMaterial('Zone4','copper');

%% Solver
s = IHLECTL3(m);clear m;
s.units.length = 'mm';
s.units.currentDensity = 'A/mm^2';
s.units.magneticVectorPotential = 'A/m';
s.setDepth(1000);
s.setGlobalFrequency(60);

s.defWinding('win','solid','voltage','series',10,0);
s.defCoil('Zone2');
s.defCoil('Zone4',1, 'negative');
s.addCoil2Winding('win', 'Zone2');
s.addCoil2Winding('win', 'Zone4');

% getting index fo boundary conditions
s.m.strefine;
s.m.strefine;
s.m.strefine;
s.m.strefine;
s.m.ggmesh;
s.clearallbcs;
s.setdbc(s.m.getfbn,0,0);

s.solve
s.evalSolidLoss
