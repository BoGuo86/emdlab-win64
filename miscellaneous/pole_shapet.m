


materialdir = [cd,'\MaterialsData'];
g = GDBC2D;

p.g = 0.5;
p.W = 30;
p.H = 40;
p.Wp = 10;
p.Hp = 32;
p.Wpsh = 8;
y = [0.5,0.6,0.8,1];


g = g.newdlinewdkps([0,p.H],[p.W,p.H],'maxLength',1);
g = g.newdlinewdkps([p.W,p.H],[p.W,0],'maxLength',1);
g = g.newdlinewdkps([p.W,0],[p.Wp+p.Wpsh,0],'maxLength',1);
g = g.newdlinewdkps([p.Wp+p.Wpsh,0],[p.Wp,0],'maxLength',1);
g = g.newdlinewdkps([p.Wp,0],[p.Wp,p.Hp],'maxLength',1);
g = g.newdlinewdkps([p.Wp,p.Hp],[p.Wp+p.Wpsh,p.Hp],'maxLength',1);
g = g.newdlinewdkps([p.Wp+p.Wpsh,0],[p.Wp+p.Wpsh,p.Hp],'maxLength',1);

y = [p.g,y];
y = p.H-y;
w = linspace(0,p.Wp+p.Wpsh,length(y));
g = g.newdlinewdkps([0,p.H-p.g],[0,p.H],'maxLength',1);

g = g.newdlinewdkps([w(1),y(1)],[w(2),y(2)],'maxLength',1);
g = g.newdlinewdkps([w(2),y(2)],[w(3),y(3)],'maxLength',1);
g = g.newdlinewdkps([w(3),y(3)],[w(4),y(4)],'maxLength',1);
g = g.newdlinewdkps([w(4),y(4)],[w(5),y(5)],'maxLength',1);

g = g.newdlinewdkps([w(5),y(5)],[w(5),p.Hp],'maxLength',1);

g = g.newcb('c','L7',1,'L6',-1,'L5',-1,'L4',-1);
g = g.newcb('a','L1',-1,'L8',-1,'L9',1,'L10',1,...
    'L11',1,'L12',1,'L13',1,'L7',-1,...
    'L3',-1,'L2',-1);

g = g.newdlinewdkps([0,p.H],[0,p.H+p.g],'maxLength',1);
g = g.newdlinewdkps([p.W,p.H],[p.W,p.H+p.g],'maxLength',1);
g = g.newdlinewdkps([0,p.H+p.g],[p.W,p.H+p.g],'maxLength',1);
g = g.newcb('g','L1',1,'L15',1,'L16',-1,'L14',-1);

g.plotbs
g = g.newdd('c',1,'c');
g = g.newdd('a',1,'a');
g = g.newdd('g',1,'g');

m = MDBCT(g);

m = m.addMaterial(materialdir,'air');

% m = m.cmirrormz('g1','g',[0,1]);
% m = m.cmirrormz('c1','c',[0,1]);
% m = m.cmirrormz('a1','a',[0,1]);
%% solver setting
s = IHNLNRMSTL3(m);clear m;
% setting units
s.scs.l = 1e-3;
s.scs.f = 1e6;
%% proccess

s = s.setExcitation('c',-1);
% s = s.setExcitation('c1',10);

s.m = s.m.ggmesh;
s.m = s.m.strefine;

% k = [s.m.getIndexOnRay([-p.W,0],[-p.W,1])
%     s.m.getIndexOnRay([p.W,0],[p.W,1])];

k0 = s.m.getIndexOnRay([0,0],[0,1]);

% [km,ks] = s.m.splitShift(k,[2*p.W,0]);
% s.m.showmeshfb;
% hold on;plot(s.m.p(k0,1),s.m.p(k0,2),'*','color','g');
% hold on;plot(s.m.p(km,1),s.m.p(km,2),'*','color','g');
% hold on;plot(s.m.p(ks,1),s.m.p(ks,2),'*','color','g');

s.m = s.m.evalKeFeC('TL3');
s = s.clearallbcs;
s = s.setdbc(unique(k0(:)),0);
% s = s.setopbc(km,ks);
s = s.assignEdata;
s = s.solve(1e-8,30);
s.plotBmagw
pause(0.1)

y = s.getbonb('g','L1');
close all
plot(y(:,2))

y = fft(y(:,2));

y = abs(y);
y = sqrt(sum(y(3:15).^2))/y(2)




