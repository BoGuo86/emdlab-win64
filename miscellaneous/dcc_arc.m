clc
clear
materialdir = [cd,'\materials'];

tmp = TMZPC([1,2,3;1,3,4],[10,0;11,0;11,1;10,1]);
tmp.strefine;
tmp.strefine;
m = TTMDBC;
m.addmz('arc', tmp.getRotateY(-90,50));
m.addMaterial(materialdir,'copper');
m.setMaterial('arc','copper');
m.ggmesh;
m.addFacetNamedSelection('in', m.getfbfiop([0,0,0],[0,0,1]));
m.addFacetNamedSelection('out', m.getfbfiop([0,0,0],[1,0,0]));
% m.shownfs

s = IHLDCCTTL4_NODAL(m);
s.setUnit('length', 'm');
s.setUnit('current', 'A');

s.setFaceCurrentBC(m.facets(m.facetNamedSelections.in,1:3),1);
tmp = m.facets(m.facetNamedSelections.out,1:3);
s.setFixedVoltageBC(unique(tmp(:)),0);

s.solve
s.plotVoltage
