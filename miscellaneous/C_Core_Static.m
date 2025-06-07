%% initialization
clear all
clc
materialdir = [cd,'\MaterialsData'];
%% Mesh
m = TMDBC;
% adding model materials
m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'m19_24ga');
m.addMaterial(materialdir,'copper');
% reading mesh file
m.read_msh_bin('C:\Users\AliJamalifard\Desktop\gmesh\c_core_s.msh');
% specifing mesh zone materials
m.setMaterial('Zone1','m19_24ga');
m.setMaterial('Zone2','copper');
m.setMaterial('Zone4','copper');
%% Getting Solver
s = IHNLNRMSTL3(m);
s.setUnit('length', 'mm');
s.setUnit('currentDensity', 'A/mm^2');
s.setUnit('magneticVectorPotential', 'A/m');
s.setDepth(1000);
%% Boundary condition
s.m.ggmesh;
s.bcs.clearAllBCs;
s.bcs.setDirichlet(s.m.getfbn,0);
%% Setting Excitations and solve
Current = 100;
Flux = zeros(1, length(Current));
for i = 1:length(Current)
    s.setExcitation('Zone2', 'value', Current(i), 'type', 'current');
    s.setExcitation('Zone4', 'value', -Current(i), 'type', 'current');
    s.solve
    Flux(i) = s.evallf('Zone2') - s.evallf('Zone4');
end
% plot(Current, Flux);
