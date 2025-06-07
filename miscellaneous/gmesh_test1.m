clc
clear


m = TTMDBC;
materialdir = [cd,'\MaterialsData'];
model = createpde;
% importGeometry(model,'Part1.stl');
importGeometry(model,'ForearmLink.stl');

% importGeometry(model,'BracketTwoHoles.stl');
% pdegplot(model,'FaceLabels','on');
generateMesh(model,'GeometricOrder','linear','Hmax',5);
m.addmz('s',TTMZPC(model.Mesh.Elements',model.Mesh.Nodes'));
m.setmzc('s', 'g');
m.ggmesh
m.setMaterial('s','copper');
m.addMaterial(materialdir,'copper');
s = NIHLTHSTTL4(m);
s.setUnit('length', 'm');
s.setUnit('volumeLossDensity', 'w/m^3');

% m.showmzs
% k0 = s.m.getfbn;
% k = s.m.getnIndexOnLine([0,0],[1,0]);
% k0 = setdiff(k0,k);
% 
% % s.setFixedTemperaturesBC(k0, 0);
% % s.setFixedTemperatureBC(k, 1);
s.setConvectionBC(s.m.facets(s.m.bfacets,1:3),.1,25);
s.setInternalLoss('s',1000000,'l');
% % s.setInternalLoss('s2',1,'l');
% % s.setInternalLoss('s3',1,'l');
% 
s.solve;

% k1 = m.getnIndexOnPlane([0,0,0],[1,0,0]);
% k2 = m.getnIndexOnPlane([0,0,60],[0,0,1]);

% m.addFacetNamedSelection('in',m.getbfiop([0,50,0],[0,0,1]));
% m.addFacetNamedSelection('out',m.getbfiop([0,0,0],[1,0,0]));
% m.shownfs
% hold on
% plot3(m.nodes(k1,1),m.nodes(k1,2),m.nodes(k1,3),'*');
% plot3(m.nodes(k2,1),m.nodes(k2,2),m.nodes(k2,3),'*');
