%% initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
% inner radius of stator
p.Rsi = 50;
% slot depth
p.hss = 10;
% ration of slot opening to slot pitch angle
p.alpha = 0.5;
% stator lamination back iron
p.bsy = 8;
% stack length of lamination
p.Lstk = 100;
% rated speed of motor 
p.rpm = 500;
% rated angular speed of motor
p.wm = p.rpm*2*pi/60;
% number of pols
p.p = 12;
% number of slots
p.Q = p.p*3;
% air gap length
p.g = 0.3;
% magnet height
p.hm = 4;
% shaft rqdius
p.rsh = 35;
% mesh size
p.msize = 0.5;
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
p8 = pmove(p5,'r',p.hm);
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
g = g.newdlinewdkps(p5,p8,p.msize);
g = g.newdarccppwdkps([0,0],p3,p8,1,p.msize);
g = g.cmirrorbd('L6',[cos(pi/p.p),sin(pi/p.p)]);
g = g.cmirrorbd('A7',[cos(pi/p.p),sin(pi/p.p)]);
g = g.newdarccpp('kp10','kp15',[0,0],1,p.msize);
g = g.newcb('rap1','L6',-1,'A9',1,'L9',-1,'A7',-1);
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
m = m.cmirrormz('c21','c11',[cos(pi/p.Q),sin(pi/p.Q)]);
m = m.setmzColor('c21','c');
for i = 1:2:2*(p.Q/p.p-1)
    m = m.crotatemz(['s',num2str(i+2)],['s',num2str(i)],2*pi/p.Q);
    m = m.crotatemz(['s',num2str(i+3)],['s',num2str(i+1)],2*pi/p.Q);
end
for i = 1:p.Q/p.p-1
    m = m.crotatemz(['c1',num2str(i+1)],['c1',num2str(i)],2*pi/p.Q);
    m = m.crotatemz(['c2',num2str(i+1)],['c2',num2str(i)],2*pi/p.Q);
end
temp = getlist('s',1:2*p.Q/p.p);
m = m.joinmzs('stator',temp{:});
m = m.cmirrormz('rap2','rap1',[1,0]);
m = m.cmirrormz('r2','r1',[1,0]);
m = m.joinmzs('rotor','r1','r2');
m = m.cmirrormz('m2','m1',[1,0]);
m = m.joinmzs('magnet','m1','m2');
m = m.ggmesh;
m.showmzs
Nts = sum(m.tzi(:,m.mzs.stator.zi));
Ntr = sum(m.tzi(:,m.mzs.rotor.zi));
%% creation of solver
s = IHNLNRMSTL3(m);clear m;
s.scs.l = 1e-3;
s.scs.f = 1e6;

simulationAngle = 2*pi/p.p;
Nsim = 50;

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
s = s.setMagnetization('magnet',@myMag1);

%% solvings
for i = 1:Ntheta
    
    s.m = s.m.ggmesh;
   
    k1 = s.m.getIndexOnCircle([0,0],p.Rsi-p.g);
    k2 = s.m.getIndexOnCircle([0,0],p.Rsi);
    s.m = s.m.makeAAG(k1,k2,p.msize);
    
     k0 = [s.m.getIndexOnCircle([0,0],p.Rsi+p.hss+p.bsy);...
        s.m.getIndexOnCircle([0,0],p.rsh)];
    k = s.m.getfb;
    k = setdiff(unique(k(:)),k0);
    [km,ks] = s.m.splitPeriodic(k,2*pi/p.p);
    
        s.m.showmeshfb
        hold on;plot(s.m.p(km,1),s.m.p(km,2),'*','color','r');
        hold on;plot(s.m.p(ks,1),s.m.p(ks,2),'*','color','c');
        hold on;plot(s.m.p(k0,1),s.m.p(k0,2),'*','color','g');
    M(i) = getframe(gcf);
    
    s = s.clearallbcs;
    s = s.setdbc(k0,0);
    s = s.setopbc(km,ks);
    s.m = s.m.evalKeFeC('TL3');
    s = s.assignEdata;
    s = s.solve(1e-6,10);
%     Torque(i) = s.evalMSTOnAG;

aLFlux(i) = s.evalLF('c12')+s.evalLF('c22');
bLFlux(i) = s.evalLF('c12')-s.evalLF('c23');
cLFlux(i) = -s.evalLF('c21')-s.evalLF('c22');

    s.plotBmag('rotor','stator','magnet')

        Bxs(:,i) = s.B(s.m.tzi(:,s.m.mzs.stator.zi),1);
        Bys(:,i) = s.B(s.m.tzi(:,s.m.mzs.stator.zi),2);
    
        Bxr(:,i) = s.B(s.m.tzi(:,s.m.mzs.rotor.zi),1);
        Byr(:,i) = s.B(s.m.tzi(:,s.m.mzs.rotor.zi),2);
    
    F(i) = getframe(gcf);
    
    if i<Ntheta
        s.m = s.m.rotatemz('rotor',xrotorAngle(2));
        s.m = s.m.rotatemz('rap1',xrotorAngle(2));
        s.m = s.m.rotatemz('rap2',xrotorAngle(2));
        s.m = s.m.rotatemz('magnet',xrotorAngle(2));
    end
    
end
close all
hold all
plot(aLFlux,'marker','o')
plot(bLFlux,'marker','d')
plot(cLFlux,'marker','s')
legend('A','B','C')
