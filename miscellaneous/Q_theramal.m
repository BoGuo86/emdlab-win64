clc
clear
materialdir = [cd,'\MaterialsData'];
[x,y] = meshgrid(linspace(0,1,10));
t = delaunay([x(:),y(:)]);

m = TMDBC();
m.addmz('s',TMZPC(t,[x(:),y(:)]));
m.addMaterial(materialdir,'copper');
m.setMaterial('s','copper');

m.ggmesh;
m.gd2elements;



s = IHLTHSTL6(m);
s.units.length = 'm';
s.units.surfaceLossDensity = 'w/m^2';

k0 = s.m.getfbn;
k = s.m.getnIndexOnLine([0,0],[1,0]);
k0 = setdiff(k0,k);

% s.m.showmd;hold on;
% plot(m.nodes(k,1),s.m.nodes(k,2),'*','color','r')
% plot(m.nodes(k0,1),s.m.nodes(k0,2),'*','color','b')

s.setdbc(k0,0);
s.setdbc(k,0);
% s.setcbc(s.m.edges(s.m.bedges,1:2),1,25);
s.setExcitation('s',1);


s.solve;
s.plotT4
