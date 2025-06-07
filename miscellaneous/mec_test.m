%% initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
%% Mesh
m = mec_dbc;
m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'iron');
m.addMaterial(materialdir,'m19_24ga');

m.addNode('n1');
m.addNode('n2');
m.addNode('ng');
m.setDatumNode('ng');
m.addReluctance('r1', 'n1', 'n2', 'material', 'iron');
m.addReluctance('r2', 'n2', 'ng', 'material', 'iron');
m.addReluctance('r3', 'n2', 'ng', 'material', 'iron', 'cArea', 1e-3);
m.addPotentialSource('is', 'ng', 'n1', 'value', 100);
clc
m.solve
m.U
