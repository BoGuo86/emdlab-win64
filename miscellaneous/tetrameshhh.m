clc
clear
materialdir = [cd,'\MaterialsData'];
h = 0.5;
xmin = -20;
xmid = -15:h:-1;

for i = 1:length(xmid)
[x,y,z] = meshgrid(xmin:h:5,1:h:2,-0.5:h:0.5);
p = [x(:),y(:),z(:)];
e = delaunay(p);
b1 = TTMZPC(p,e);

[x,y,z] = meshgrid(xmin:h:5,-1:-h:-2,-0.5:h:0.5);
p = [x(:),y(:),z(:)];
e = delaunay(p);
b2 = TTMZPC(p,e);

[x,y,z] = meshgrid(xmid(i):h:xmid(i)+1,-3:h:3,0.5:h:1.5);
p = [x(:),y(:),z(:)];
e = delaunay(p);
b3 = TTMZPC(p,e);

m = TTMDBC('b1',b1,'b2',b2,'b3',b3);

m = m.addMaterial(materialdir,'air');

s = IHLTHSTTL4(m);

s.m = s.m.ggmesh;


k = s.m.getIndexOnPlane([xmin,0,0],[1,0,0]);

[k1,k2] = s.m.splitShift(k,[0,3,0]);


s = s.setdbc(k1,1);
s = s.setdbc(k2,0);

% s = s.setExcitation('b3',0.025);

s.m = s.m.evalKeFeC('TTL4');
s = s.assignEdata;

s = s.solve;
% s = s.evalQ;
s.plotT

% set(gcf,'units','normalized','outerposition',[0 0 1 1])
F(i) = getframe(gcf);

end





