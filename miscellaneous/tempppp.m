

glib_test3;
close all
m = TTMDBC;
s.m.getMeshZoneExtrude(m,'stator',0:2:14);
for i = 1:36
s.m.getMeshZoneExtrude(m,['c1',num2str(i)],0:2:20);
s.m.getMeshZoneExtrude(m,['c2',num2str(i)],0:2:20);
end
% s.m.getMeshZoneExtrude(m,'rotor',-20:2:20);
% m.showmzs
m.ggmesh;
m.addFacetNamedSelection('Rear', m.getbfiop([0,0,0],[0,0,1]));
m.addFacetNamedSelection('Front', m.getbfiop([0,0,20],[0,0,1]));
m.shownfs
