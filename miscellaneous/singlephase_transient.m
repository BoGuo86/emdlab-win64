%% initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
g = GDBC2D;
% parameters
p.wc = 5;
p.wy = 5;
p.hy = 20;
p.ly = 20;
p.g = 0.5;
p.msy = 1.5;
p.msgin = 1.5;
p.msgout = 1.5;
%% creation of yoke
yp = [0,0
    p.wy,0
    p.wy,p.hy-p.wy
    p.wy+p.wc,p.hy-p.wy
    p.ly,p.hy-p.wy
    p.ly,p.hy
    0,p.hy];
g = g.newpolygonald('yoke',p.msy,yp,p.msy);
g = g.newdlinewdkps([p.wy,0],[p.wy+p.wc,0],p.msgin);
%% creation of coil
g = g.newdlinewdkps([p.wy+p.wc,0],[p.wy+p.wc,p.hy-p.wy],p.msgin);
g = g.newcbwd('coil',p.msgin,'L8',1,'L9',1,'L3',-1,'L2',-1);
%% creation of armature
ap = [p.ly+p.g,0
    p.ly+p.g+p.wy,0
    p.ly+p.g+p.wy,p.hy
    p.ly+p.g,p.hy];
g = g.newpolygonald('armature',p.msy,ap,p.msy);
%% creation of air1 and air 2
g = g.newdlinewdkps([p.wy+p.wc,0],[p.ly+p.g,0],p.msgin);
g = g.newdlinewdkps([p.wy+p.ly+p.g,0],[3*p.ly,0],p.msgout);
g = g.newdlinewdkps([3*p.ly,0],[3*p.ly,3*p.hy],p.msgout);
g = g.newdlinewdkps([3*p.ly,3*p.hy],[0,3*p.hy],p.msgout);
g = g.newdlinewdkps([0,3*p.hy],[0,p.hy],p.msgout);
g = g.newdlinewdkps([p.ly,p.hy],[p.ly+p.g,p.hy],p.msgin);
g = g.newcbwd('airout',p.msgout,'L15',1,'L16',1,'L17',1,'L18',1,...
    'L6',-1,'L19',1,'L12',-1,'L11',-1);
g = g.newcbwd('airin',p.msgin,'L14',1,'L13',-1,'L19',-1,'L5',-1,...
    'L4',-1,'L9',-1);
%% generation of mesh
m = MDBCT(g);clear g;
m = m.setMaterial('yoke','m19_24ga');
m = m.setMaterial('armature','m19_24ga');
m = m.setMaterial('coil','copper');
m = m.addMaterial(materialdir,'air');
m = m.addMaterial(materialdir,'m19_24ga');
m = m.addMaterial(materialdir,'copper');
%% setting solve
s = IHNLNRMDCFWTMTL3(m);clear m;
s.m = s.m.ggmesh;
k0 = [s.m.getindex('airout','L16')
    s.m.getindex('airout','L17')
    s.m.getindex('airout','L18')
    s.m.getindex('yoke','L7')];
s.m = s.m.evalKeFeC('TL3');
s.m = s.m.evalMe;
s.scs.l = 1e-3;
s.scs.f = 1e6;
s.scs.t = 1;
s = s.setdbc(k0,0);
%% solvings
ts = 1e-3;
tf = 100e-3;
Nt = tf/ts;
t = 0:ts:20e-3;
s = s.setExcitation('coil',t,sin(2*pi*50*t),100);
s = s.solve(1e-6,10,0.001,ts,Nt);
plot(0:ts:ts*Nt,[0,s.linkageFluxes(2,:)]*100)
% hold on;plot(0:ts:ts*Nt,[0,s.currents(2,:)]/max([0,s.currents(2,:)]))





