%% Steady State Simulation And Calculation of Phase Currents Versus Postion 
% on and off angle (mechanical degree)
clc
clear
rpm = 1500;
f = fopen('C:\Users\AliJamalifard\Desktop\ia.tab','r');
i_phaseA = fscanf(f,'%f');
fclose(f);
i_phaseA = reshape(i_phaseA,2,[]);
i_phaseA(1,:) = rpm*(pi/30)*i_phaseA(1,:)*1e-3;
simulationAngle = pi/2;
Nsim = 1;
%% initialization
clc
materialdir = [cd,'\MaterialsData'];
g = GDBC2D;
% parameters
p.rsh = 10.000000;
p.R0 = 30.400000;
p.R1 = 40.900000;
p.R2 = 65.800000;
p.R3 = 75.00000;
p.betar = 19 * pi/180;
p.betas = 16 * pi/180;
p.g = 0.3;
p.dw = 6.8;
p.Rw = p.R1+p.g+p.dw;
p.gammar = asin(p.R1*sin(p.betar/2)/p.R0);
p.gammas = asin((p.R1+p.g)*sin(p.betas/2)/p.R2);
p.gammaw = asin((p.R1+p.g)*sin(p.betas/2)/p.Rw);
% mesh size
p.hr = 1;
p.hs =1;
p.hg =1;
p.thetag = 1.5;
%% creation of rotor
p1 = [p.rsh,0];
p2 = [p.R1,0];
p3 = p.R1*[cos(p.betar/2),sin(p.betar/2)];
p4 = p.R0*[cos(p.gammar),sin(p.gammar)];
p5 = p.R0*[cos(pi/8),sin(pi/8)];
p6 = p.rsh*[cos(pi/8),sin(pi/8)];
g = g.newdlinewdkps(p1,p2,'maxLength',0.9*p.hr);
g = g.newdarccppwdkps([0,0],p2,p3,'maxLength',0.9*min(p.hr,p.hg),'maxDegree',p.thetag);
g = g.newdlinewdkps(p3,p4,'maxLength',0.9*min(p.hr,p.hg));
g = g.newdarccppwdkps([0,0],p4,p5,'maxLength',0.9*min(p.hr,p.hg),'maxDegree',5);
g = g.newdlinewdkps(p5,p6,'maxLength',0.9*p.hr);
g = g.newdarccppwdkps([0,0],p6,p1,'direction',-1,'maxLength',0.9*p.hr,'maxDegree',15);
g = g.newcb('r1','L1',1,'A1',1,'L2',1,'A2',1,'L3',1,'A3',1);

%% creation of stator
p1 = [p.R1+p.g,0];
p2 = [p.R3,0];
p3 = p.R3*[cos(pi/12),sin(pi/12)];
p4 = p.R2*[cos(pi/12),sin(pi/12)];
p5 = p.R2*[cos(p.gammas),sin(p.gammas)];
p6 = (p.R1+p.g)*[cos(p.betas/2),sin(p.betas/2)];
g = g.newdlinewdkps(p1,p2,'maxLength',0.9*p.hs);
g = g.newdarccppwdkps([0,0],p2,p3,'maxLength',0.9*p.hs,'maxDegree',5);
g = g.newdlinewdkps(p3,p4,'maxLength',0.9*p.hs);
g = g.newdarccppwdkps([0,0],p4,p5,'direction',-1,'maxLength',0.9*min(p.hs,p.hg),'maxDegree',5);
g = g.newdlinewdkps(p5,p6,'maxLength',0.9*min(p.hs,p.hg));
g = g.newdarccppwdkps([0,0],p6,p1,'direction',-1,'maxLength',0.9*min(p.hs,p.hg),'maxDegree',p.thetag);
g = g.newcb('s1','L4',1,'A4',1,'L5',1,'A5',1,'L6',1,'A6',1);

%% creation of coil
pw = (p.R1+p.g)*[cos(pi/12),sin(pi/12)];
g = g.newdlinewdkps(p4,pw,'maxLength',0.9*p.hg);
g = g.newdarccppwdkps([0,0],pw,p6,'direction',-1,'maxLength',0.9*p.hg,'maxDegree',p.thetag);
g = g.newcb('c1','L6',-1,'A5',-1,'L7',1,'A7',1);

%% creation of rotor air pockets
g = g.cmirrorbd('L2',[cos(pi/8),sin(pi/8)]);
g = g.cmirrorbd('A2',[cos(pi/8),sin(pi/8)]);
g = g.newdarccpp('kp3','kp14',[0,0],'maxLength',0.95*p.hg,'maxDegree',p.thetag);
g = g.newcb('rap1','L2',-1,'A9',1,'L8',1,'A8',1,'A2',-1);

% g = g.newdd('r1',p.hr,'r1');
% g = g.newdd('s1',p.hs,'s1');
% g = g.newdd('c1',p.hg,'c1');
% g = g.newdd('rap1',p.hg,'rap1');

g = g.newdDM('r1','r1');
g = g.newdDM('s1','s1');
g = g.newdDM('c1','c1');
g = g.newdDM('rap1','rap1');

%% creation of mesh
m = MDBCT(g);clear g;
m = m.setmzColor('r1','r');
m = m.setmzColor('s1','r');
m = m.setmzColor('c1','m');
m = m.setmzColor('rap1','w');
m = m.setMaterial('r1','m19_24ga');
m = m.setMaterial('s1','m19_24ga');
m = m.setMaterial('c1','copper');
m = m.addMaterial(materialdir,'air');
m = m.addMaterial(materialdir,'m19_24ga');
m = m.addMaterial(materialdir,'copper');
m = m.cmirrormz('r2','r1',[cos(pi/8),sin(pi/8)]);
m = m.crotatemz('r3','r1',pi/4);
m = m.crotatemz('r4','r2',pi/4);
m = m.cmirrormz('s2','s1',[cos(pi/12),sin(pi/12)]);
m = m.crotatemz('s3','s1',pi/6);
m = m.crotatemz('s4','s2',pi/6);
m = m.crotatemz('s5','s3',pi/6);
m = m.crotatemz('s6','s4',pi/6);
m = m.cmirrormz('c2','c1',[cos(pi/12),sin(pi/12)]);
m = m.crotatemz('c3','c1',pi/6);
m = m.crotatemz('c4','c2',pi/6);
m = m.crotatemz('c5','c3',pi/6);
m = m.crotatemz('c6','c4',pi/6);
m = m.crotatemz('rap2','rap1',pi/4);
%% joining mesh
temp = getlist('r',1:4);
m = m.joinmzs('rotor',temp{:});
temp = getlist('s',1:6);
m = m.joinmzs('stator',temp{:});
m = m.rotatemz('rotor',-pi/8);
m = m.rotatemz('rap1',-pi/8);
m = m.rotatemz('rap2',-pi/8);
m = m.ggmesh;
szi = m.mzs.stator.zi;
rzi = m.mzs.rotor.zi;
Nts = sum(m.tzi(:,szi));
Ntr = sum(m.tzi(:,rzi));
%% solver setting
s = IHNLNRMSTL3(m);clear m;
s.scs.l = 1e-3;
s.scs.f = 1e6;
%% proccess
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

temp = 0;
for i = 1:Ntheta
    
    s.m = s.m.ggmesh;
    kin = s.m.getIndexOnCircle([0,0],p.R1);
    kout = s.m.getIndexOnCircle([0,0],p.R1+p.g);
    s.m = s.m.makeAAG(kin,kout,1,p.thetag);
    % getting index fo boundary conditions
    k0 = [s.m.getIndexOnCircle([0,0],p.rsh);...
        s.m.getIndexOnCircle([0,0],p.R3)];
    k = s.m.getfb;
    k = setdiff(unique(k(:)),k0);
    [km,ks] = s.m.splitPeriodic(k,pi/2);
    
%         s.m.showmeshfb
%         hold on;plot(s.m.p(km,1),s.m.p(km,2),'*','color','r');
%         hold on;plot(s.m.p(ks,1),s.m.p(ks,2),'*','color','c');
%         hold on;plot(s.m.p(k0,1),s.m.p(k0,2),'*','color','g');
%     M(i) = getframe(gcf);
    
    s.m = s.m.evalKeFeC('TL3');
    s = s.clearallbcs;
    s = s.setdbc(k0,0);
    s = s.setopbc(km,ks);
    
    ia(i) = interp1(i_phaseA(1,:),i_phaseA(2,:),rem(xrotorAngle(i),pi/4));
    ib(i) = interp1(i_phaseA(1,:),i_phaseA(2,:),rem(xrotorAngle(i)+15*pi/180,pi/4));
    ic(i) = interp1(i_phaseA(1,:),i_phaseA(2,:),rem(xrotorAngle(i)+30*pi/180,pi/4));
    
    s = s.setExcitation('c1',ia(i)*188,'C');
    s = s.setExcitation('c6',ia(i)*188,'C');
    s = s.setExcitation('c2',ib(i)*188,'C');
    s = s.setExcitation('c3',-ib(i)*188,'C');
    s = s.setExcitation('c4',-ic(i)*188,'C');
    s = s.setExcitation('c5',ic(i)*188,'C');
    
    s = s.assignEdata;
    s = s.solve(1e-6,10);
    s = s.evalB;
    
    aLFlux(i) = s.evalLF('c1')+s.evalLF('c6');
    bLFlux(i) = s.evalLF('c2')-s.evalLF('c3');
    cLFlux(i) = s.evalLF('c5')-s.evalLF('c4');
    
    
%     s.plotBmagSmooth
%     hold on
%     s.m.plotwf
    Bxs(:,i) = s.B(s.m.tzi(:,szi),1);
    Bys(:,i) = s.B(s.m.tzi(:,szi),2);
    
    Bxr(:,i) = s.B(s.m.tzi(:,rzi),1);
    Byr(:,i) = s.B(s.m.tzi(:,rzi),2);
    
%     F(i) = getframe(gcf);
    
    if i<Ntheta
        s.m = s.m.rotatemz('rotor',xrotorAngle(2));
        s.m = s.m.rotatemz('rap1',xrotorAngle(2));
        s.m = s.m.rotatemz('rap2',xrotorAngle(2));
    end

    
end

aLFlux = aLFlux*2*188*55/1000;
bLFlux = bLFlux*2*188*55/1000;
cLFlux = cLFlux*2*188*55/1000;

close all
subplot(211)
hold all
plot(xrotorAngle*180/pi,ia)
plot(xrotorAngle*180/pi,ib)
plot(xrotorAngle*180/pi,ic)
xlabel('Rotor Angle [Deg]')
ylabel('Phase Current [A]')
legend('Phase A','Phase B','Phase C')
subplot(212)
hold all
plot(xrotorAngle*180/pi,aLFlux)
plot(xrotorAngle*180/pi,bLFlux)
plot(xrotorAngle*180/pi,cLFlux)
xlabel('Rotor Angle [Deg]')
ylabel('Phase Linkage Flux [wb]')
legend('Phase A','Phase B','Phase C')