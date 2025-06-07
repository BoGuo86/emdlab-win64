%% initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
% inner radius of stator
p.Rsi = 50;
% slot depth
p.hss = 10;
% ration of slot opening to slot pitch angle
p.alpha = 0.6;
% stator lamination back iron
p.bsy = 8;
% number of slots
p.Q = 48;
% number of pols
p.p = 12;
% air gap length
p.g = 0.4;
% magnet height
p.hm = 4;
% shaft rqdius
p.rsh = 35;
p.msize = 2;
%% creation of geometry
g = GDBC2D;
% stator lamination
p1 = [p.Rsi,0];
p2 = [p.Rsi+p.hss+p.bsy,0];
p3 = pmove(p2,'theta',pi/p.Q);
p4 = pmove(p3,'r',-p.bsy);
p5 = pmove(p4,'theta',-p.alpha*pi/p.Q);
p6 = pmove(p5,'r',-p.hss);
p7 = pmove(p4,'r',-p.hss);
g = g.newdlinewdkps(p1,p2,p.msize);
g = g.newdarccppwdkps([0,0],p2,p3,1,p.msize);
g = g.newdlinewdkps(p3,p4,p.msize);
g = g.newdarccppwdkps([0,0],p4,p5,-1,p.msize);
g = g.newdlinewdkps(p5,p6,p.msize);
g = g.newdarccppwdkps([0,0],p6,p1,-1,p.msize);
g = g.newcb('s1','L1',1,'A1',1,'L2',1,'A2',1,'L3',1,'A3',1);
% coils
g = g.newdlinewdkps(p4,p7,p.msize);
g = g.newdarccppwdkps([0,0],p6,p7,1,p.msize);
g = g.newcb('c11','L3',-1,'A2',-1,'L4',1,'A4',-1);
% creation of magnet
p1 = [p.Rsi-p.g-p.hm,0];
p2 = [p.Rsi-p.g,0];
p3 = pmove(p2,'theta',2*pi/p.p/3);
p4 = pmove(p3,'r',-p.hm);
p5 = pmove(p4,'theta',1*pi/p.p/3);
p7 = [p.rsh,0];
p6 = pmove(p7,'theta',pi/p.p);
g = g.newdlinewdkps(p1,p2,p.msize);
g = g.newdarccppwdkps([0,0],p2,p3,1,p.msize);
g = g.newdlinewdkps(p3,p4,p.msize);
g = g.newdarccppwdkps([0,0],p4,p1,-1,p.msize);
g = g.newcb('m1','L5',1,'A5',1,'L6',1,'A6',1);
% rotor lamination
g = g.newdarccppwdkps([0,0],p4,p5,1,p.msize);
g = g.newdlinewdkps(p5,p6,p.msize);
g = g.newdarccppwdkps([0,0],p6,p7,-1,p.msize,5);
g = g.newdlinewdkps(p7,p1,p.msize);
g = g.newcb('r1','L8',1,'A6',-1,'A7',1,'L7',1,'A8',1);
% rotor air pocket
g = g.cmirrorbd('L6',[cos(pi/p.p),sin(pi/p.p)]);
g = g.cmirrorbd('A7',[cos(pi/p.p),sin(pi/p.p)]);
g = g.newdarccpp('kp10','kp15',[0,0],1,p.msize);
g = g.newcb('rap1','L6',-1,'A10',1,'L9',1,'A9',1,'A7',-1);
% creation of domains
g = g.newdDM('s1',p.msize,'s1');
g = g.newdDM('c11',p.msize,'c11');
g = g.newdDM('m1',p.msize,'m1');
g = g.newdDM('r1',p.msize,'r1');
g = g.newdDM('rap1',p.msize,'rap1');
close all
%% creation of mesh
m = MDBCT(g);clear g
m = m.setMaterial('r1','m19_24ga');
m = m.setMaterial('s1','m19_24ga');
m = m.setMaterial('c11','copper');
m = m.addMaterial(materialdir,'air');
m = m.addMaterial(materialdir,'m19_24ga');
m = m.addMaterial(materialdir,'copper');
m = m.setmzColor('s1','r');
m = m.setmzColor('r1','r');
m = m.setmzColor('c11','g');
m = m.setmzColor('m1','m');
m = m.setmzColor('rap1','w');
m = m.cmirrormz('s2','s1',[cos(pi/p.Q),sin(pi/p.Q)]);
m = m.cmirrormz('r2','r1',[1,0]);
m = m.cmirrormz('c21','c11',[cos(pi/p.Q),sin(pi/p.Q)]);
m = m.setmzColor('c21','c');
for i = 1:2:2*(p.Q-1)
    m = m.crotatemz(['s',num2str(i+2)],['s',num2str(i)],2*pi/p.Q);
    m = m.crotatemz(['s',num2str(i+3)],['s',num2str(i+1)],2*pi/p.Q);
end
for i = 1:2:2*(p.p-1)
    m = m.crotatemz(['r',num2str(i+2)],['r',num2str(i)],2*pi/p.p);
    m = m.crotatemz(['r',num2str(i+3)],['r',num2str(i+1)],2*pi/p.p);
end
for i = 1:p.Q-1
    m = m.crotatemz(['c1',num2str(i+1)],['c1',num2str(i)],2*pi/p.Q);
    m = m.crotatemz(['c2',num2str(i+1)],['c2',num2str(i)],2*pi/p.Q);
end
temp = getlist('s',1:2*p.Q);
m = m.joinmzs('stator',temp{:});
temp = getlist('r',1:2*p.p);
m = m.joinmzs('rotor',temp{:});
m = m.cmirrormz('m2','m1',[1,0]);
m = m.joinmzs('magnet1','m1','m2');
for i = 1:p.p-1
    m = m.crotatemz(['magnet',num2str(i+1)],['magnet',num2str(i)],2*pi/p.p);
    m = m.crotatemz(['rap',num2str(i+1)],['rap',num2str(i)],2*pi/p.p);
end
m = m.ggmesh;
Nts = sum(m.tzi(:,74));
Ntr = sum(m.tzi(:,75));
%% creation of solver
s = IHNLNRMSTL3(m);clear m;
s.scs.l = 1e-3;
s.scs.f = 1e6;
%%
simulationAngle = 2*2*pi/p.p;
Nsim = 150;

xrotorAngle = linspace(0,simulationAngle,Nsim);
Ntheta = length(xrotorAngle);
F(Ntheta) = struct('cdata',[],'colormap',[]);
M(Ntheta) = struct('cdata',[],'colormap',[]);

aLFlux = zeros(Ntheta,1);
bLFlux = zeros(Ntheta,1);
cLFlux = zeros(Ntheta,1);

Torque = zeros(Ntheta,1);
Bxs = zeros(Nts,Ntheta);
Bys = zeros(Nts,Ntheta);
Bxr = zeros(Ntr,Ntheta);
Byr = zeros(Ntr,Ntheta);

ia = zeros(Ntheta,1);
ib = zeros(Ntheta,1);
ic = zeros(Ntheta,1);

% setting magnetization of magnets
for i = 1:p.p
    if rem(i,2) == 0
        s = s.setMagnetization(['magnet',num2str(i)],@myMag1);
    else
        s = s.setMagnetization(['magnet',num2str(i)],@myMag2);
    end
end

%% solvings
for i = 1:Ntheta
    
    s.m = s.m.ggmesh;
    k0 = [s.m.getIndexOnCircle([0,0],p.Rsi+p.hss+p.bsy);...
        s.m.getIndexOnCircle([0,0],p.rsh)];
    k1 = s.m.getIndexOnCircle([0,0],p.Rsi-p.g);
    k2 = s.m.getIndexOnCircle([0,0],p.Rsi);
    s.m = s.m.makeAG(k1,k2);
    
    s = s.clearallbcs;
    s = s.setdbc(k0,0);
    s.m = s.m.evalKeFeC('TL3');
    s = s.assignEdata;
    s = s.solve(1e-6,10);
    Torque(i) = s.evalMSTOnAG;
    a = getlist('magnet',1:12);
    s.plotBmag('rotor','stator',a{:})
    
        Bxs(:,i) = s.B(s.m.tzi(:,74),1);
        Bys(:,i) = s.B(s.m.tzi(:,74),2);
    
        Bxr(:,i) = s.B(s.m.tzi(:,75),1);
        Byr(:,i) = s.B(s.m.tzi(:,75),2);
    
    F(i) = getframe(gcf);
    
    if i<Ntheta
        s.m = s.m.rotatemz('rotor',xrotorAngle(2));
        for j = 1:p.p
            s.m = s.m.rotatemz(['rap',num2str(j)],xrotorAngle(2));
            s.m = s.m.rotatemz(['magnet',num2str(j)],xrotorAngle(2));
        end
    end
    
    
    
end