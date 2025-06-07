%% initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
g = GDBC2D;
% parameters
p.Rs1 = 40;
p.Rs2 = 55;
p.bsy = 5;
p.Nss = 48;
p.hst = 10;
p.ts = 0.5;
p.hs = 1;
%% creation of stator
p1 = [p.Rs1,0];
p2 = [p.Rs2,0];
p3 = protate(p2,pi/p.Nss);
p4 = protate([p.Rs1+p.hst,0],pi/p.Nss);
p5 = protate(p1,pi*p.ts/p.Nss);
g = g.newdlinewdkps(p1,p2,p.hs);
g = g.newdarccppwdkps([0,0],p2,p3,1,p.hs);
g = g.newdarccppwdkps([0,0],p1,p5,1,p.hs);
g = g.newdlinewdkps(p3,p4,p.hs);
g = g.newdlinewdkps(p5,p5+[p.hst,0],p.hs);
g = g.newdlinewdkps(p5+[p.hst,0],p4,p.hs);
g = g.filletll('L3','L4',1,p.hs,30);
g = g.newcbwd('s1',p.hs,'L1',1,'A1',1,'L2',1,'L4',-1 ...
,'A3',-1,'L3',-1,'A2',-1);
%% creation of coil
g = g.cmirrorbd('L3',[cos(pi/p.Nss),sin(pi/p.Nss)]);
g = g.cmirrorbd('A3',[cos(pi/p.Nss),sin(pi/p.Nss)]);
g = g.cmirrorbd('L4',[cos(pi/p.Nss),sin(pi/p.Nss)]);
g = g.newdarccpp('kp4','kp9',[0,0],1,p.hs);
g = g.newcbwd('c1',p.hs,'L3',1,'A3',1,'L4',1,'L6',-1 ...
,'A4',-1,'L5',-1,'A5',-1);
%% creation of airgap
g = g.cmirrorbd('A2',[cos(pi/p.Nss),sin(pi/p.Nss)]);
g = g.newdlinewdkps([p.Rs1-3,0],[p.Rs1,0],p.hs);
g = g.newdlinewdkps(protate([p.Rs1-3,0],2*pi/p.Nss),protate([p.Rs1,0],2*pi/p.Nss),p.hs);
g = g.newdarccpp('kp13','kp14',[0,0],1,p.hs);
g = g.newcbwd('a1',p.hs,'L7',1,'A2',1,'A5',1,'A6',-1 ...
,'L8',-1,'A7',-1);
%%
m = MDBCT(g);clear g;
m = m.setMaterial('s1','m19_24ga');
% m = m.setMaterial('c1','copper');
m = m.addMaterial(materialdir,'air');
m = m.addMaterial(materialdir,'m19_24ga');
m = m.cmirrormz('s2','s1',[cos(pi/p.Nss),sin(pi/p.Nss)]);
for i = 1:2:2*(p.Nss-1)
    m = m.crotatemz(['s',num2str(i+2)],['s',num2str(i)],2*pi/p.Nss);
    m = m.crotatemz(['s',num2str(i+3)],['s',num2str(i+1)],2*pi/p.Nss); 
end
for i = 2:p.Nss
    m = m.crotatemz(['c',num2str(i)],['c',num2str(i-1)],2*pi/p.Nss);
    m = m.crotatemz(['a',num2str(i)],['a',num2str(i-1)],2*pi/p.Nss);
end
m = m.ggmesh;
m.showmzs
s = IHNLNRMDCFWTMTL3(m);clear m;
%% proccess
s.m = s.m.ggmesh;
k0 = s.m.getIndexOnCircle([0,0],p.Rs2);
% hold on;plot(s.m.p(km,1),s.m.p(km,2),'*','color','r');
% hold on;plot(s.m.p(ks,1),s.m.p(ks,2),'*','color','c');
% hold on;plot(s.m.p(k0,1),s.m.p(k0,2),'*','color','g');
s.m = s.m.evalKeFeC('TL3');
s.m = s.m.evalMe;
s.scs.l = 1e-3;
s.scs.f = 1e6;
s.scs.t = 1e-3;
s = s.clearallbcs;
s = s.setdbc(k0,0);
% s = s.setopbc(km,ks);
%% solvings
ts = 0.25;
tf = 20;
Nt = tf/ts;
t = 0:ts:20;
Np = 20;
ia = 5*sin(2*pi*50*t*1e-3);
ib = 5*sin(2*pi*50*t*1e-3+2*pi/3);
ic = 5*sin(2*pi*50*t*1e-3-2*pi/3);
for i = [1 2 3 4 25 26 27 28]
    s = s.setExcitation(['c',num2str(i)],t,ia,Np,'CF');
end
for i = [1 2 3 4 25 26 27 28] + 4
    s = s.setExcitation(['c',num2str(i)],t,-ib,Np,'CF');
end
for i = [1 2 3 4 25 26 27 28] + 8
    s = s.setExcitation(['c',num2str(i)],t,ic,Np,'CF');
end
for i = [1 2 3 4 25 26 27 28] + 12
    s = s.setExcitation(['c',num2str(i)],t,-ia,Np,'CF');
end
for i = [1 2 3 4 25 26 27 28] + 16
    s = s.setExcitation(['c',num2str(i)],t,ib,Np,'CF');
end
for i = [1 2 3 4 25 26 27 28] + 20
    s = s.setExcitation(['c',num2str(i)],t,-ic,Np,'CF');
end

s = s.solve(1e-6,10,0.001,ts,Nt);


% 
% 

