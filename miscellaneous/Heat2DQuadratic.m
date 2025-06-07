%% initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
g = GDBC2D;
% parameters
param.rsh = 10.000000;
param.R0 = 30.400000;
param.R1 = 40.900000;
param.R2 = 65.800000;
param.R3 = 75.00000;
param.betar = 19 * pi/180;
param.betas = 16 * pi/180;
param.wh = 4;
param.wc = 2;
param.g = 0.3;
param.Np = 188;
param.Rw = param.R1 + param.g + param.wc;
param.gammaw = asin((param.R1+param.g)*sin(param.betas/2)/(param.R1+param.g+param.wc));
param.gammar = asin(param.R1*sin(param.betar/2)/param.R0);
param.gammas = asin((param.R1+param.g)*sin(param.betas/2)/param.R2);
% mesh size in rotor, stator and air gap
param.hr = 2.5;
param.hs =1.5;
param.hg =1.5;
param.thetag = 3;
%% creation of rotor
p1 = [param.rsh,0];
p2 = [param.R1,0];
p3 = param.R1*[cos(param.betar/2),sin(param.betar/2)];
p4 = param.R0*[cos(param.gammar),sin(param.gammar)];
p5 = param.R0*[cos(pi/8),sin(pi/8)];
p6 = param.rsh*[cos(pi/8),sin(pi/8)];
g = g.newdlinewdkps(p1,p2,'maxLength',0.9*param.hr);
g = g.newdarccppwdkps([0,0],p2,p3,'maxLength',0.9*min(param.hr,param.hg),'maxDegree',param.thetag);
g = g.newdlinewdkps(p3,p4,'maxLength',0.9*min(param.hr,param.hg));
g = g.newdarccppwdkps([0,0],p4,p5,'maxLength',0.9*min(param.hr,param.hg),'maxDegree',5);
g = g.newdlinewdkps(p5,p6,'maxLength',0.9*param.hr);
g = g.newdarccppwdkps([0,0],p6,p1,'direction',-1,'maxLength',0.9*param.hr,'maxDegree',15);
g = g.filletla('L2','A2',2.5,'Nnodes',5);
g = g.newcb('r1','L1',1,'A1',1,'L2',1,'A4',-1,'A2',1,'L3',1,'A3',1);
%% creation of stator
p1 = [param.R1+param.g,0];
p2 = [param.R3,0];
p3 = param.R3*[cos(pi/12),sin(pi/12)];
p4 = param.R2*[cos(pi/12),sin(pi/12)];
p5 = param.R2*[cos(param.gammas),sin(param.gammas)];
p6 = (param.R1+param.g)*[cos(param.betas/2),sin(param.betas/2)];
p7 = param.Rw*[cos(param.gammaw),sin(param.gammaw)];
g = g.newdlinewdkps(p1,p2,'maxLength',0.9*param.hs);
g = g.newdarccppwdkps([0,0],p2,p3,'maxLength',0.9*param.hs,'maxDegree',5);
g = g.newdlinewdkps(p3,p4,'maxLength',0.9*param.hs);
g = g.newdarccppwdkps([0,0],p4,p5,'direction',-1,'maxLength',0.9*min(param.hs,param.hg),'maxDegree',5);
g = g.newdlinewdkps(p5,p7,'maxLength',0.9*min(param.hs,param.hg));
g = g.newdlinewdkps(p7,p6,'maxLength',0.9*min(param.hs,param.hg));
g = g.newdarccppwdkps([0,0],p6,p1,'direction',-1,'maxLength',0.9*min(param.hs,param.hg),'maxDegree',param.thetag);
g = g.newcb('s1','L4',1,'A5',1,'L5',1,'A6',1,'A8',-1,'L6',1,'L7',1,'A7',1);
g = g.filletla('L6','A6',2.5,'Nnodes',5);
%% creation of housing
p1 = [param.R3,0];
p2 = p1 + [param.wh,0];
p3 = pmove(p2,'theta',pi/12);
p4 = pmove(p1,'theta',pi/12);
g = g.newdlinewdkps(p1,p2,'maxLength',0.9*param.hs);
g = g.newdarccppwdkps([0,0],p2,p3,'maxLength',0.9*param.hs,'maxDegree',5);
g = g.newdlinewdkps(p3,p4,'maxLength',0.9*param.hs);
g = g.newcb('h1','L8',1,'A9',1,'L9',1,'A5',-1);

%% creation of coil
p1 = param.Rw*[cos(param.gammaw),sin(param.gammaw)];
p2 = pmove([param.Rw,0],'theta',pi/12);
p3 = pmove(p2,'r',param.R2-param.Rw);
g = g.newdlinewdkps(p2,p3,'maxLength',0.9*param.hg);
g = g.newdarccppwdkps([0,0],p1,p2,'maxLength',0.9*param.hg,'maxDegree',param.thetag);

g = g.newcb('c1','L6',-1,'A8',1,'A6',-1,'L10',-1,'A10',-1);

%% creation of end winding
g = g.newdlinewdkps([param.Rw,0],[param.R2,0],'maxLength',0.9*param.hg);
g = g.newdarccpp('kp21','kp14',[0,0],'maxLength',0.9*param.hg,'maxDegree',5);
g = g.newdarccpp('kp22','kp17',[0,0],'maxLength',0.9*param.hg,'maxDegree',5);
g = g.newcb('ew1','L11',1,'A12',1,'A8',-1,'L6',1,'A11',-1);
%% creation of wedge
p1 = param.Rw*[cos(pi/12),sin(pi/12)];
p2 = (param.R1+param.g)*[cos(pi/12),sin(pi/12)];
g = g.newdlinewdkps(p1,p2,'maxLength',0.9*param.hg);
g = g.newdarccpp('kp15','kp23',[0,0],'maxLength',0.9*param.hg,'maxDegree',param.thetag);
g = g.newcb('w','L7',-1,'A10',1,'L12',1,'A13',-1);

% g.plotbs

% g = g.newdDM('r1','r1');
g = g.newdDM('s1','s1');
g = g.newdDM('h1','h1');
g = g.newdDM('c1','c1');
g = g.newdDM('w','w');
%%

m = MDBCT(g);

m = m.addMaterial(materialdir,'lamination');
m = m.addMaterial(materialdir,'aluminium');
m = m.addMaterial(materialdir,'win');
m = m.addMaterial(materialdir,'wedge');

m = m.setMaterial('s1','lamination');
m = m.setMaterial('h1','aluminium');
m = m.setMaterial('c1','win');
m = m.setMaterial('w','wedge');

m = m.setmzColor('h1',[0.7617,0.7617,0.7617]);
m = m.setmzColor('c1',[0.9688,0.6016,0.2383]);
m = m.setmzColor('s1',[0.4180,0.4180, 0.7109]);
m = m.setmzColor('w',[0.3,0, 0]);

close all

%%

s = IHLTHSTL3(m);
s.scs.l = 1e-3;
s.scs.f = 1e6;
s.m = s.m.ggmesh;
s.m = s.m.strefine;
s.m = s.m.gd2elements;
s.m = s.m.setdata;
s.m.showmeshfb

kins1 = s.m.getIndexOnRay([0,0],[cos(pi/12),sin(pi/12)]);
kins2 = s.m.getIndexOnRay([0,0],[1,0]);
eout = s.m.getEindex('h1','A9');
eini = s.m.getEindex('s1','A7');
einc = s.m.getEindex('w','A13');

hold on;plot(s.m.p(kins1(:),1),s.m.p(kins1(:),2),'*')
hold on;plot(s.m.p(kins2(:),1),s.m.p(kins2(:),2),'*')
hold on;plot(s.m.p(einc(:),1),s.m.p(einc(:),2),'*')
hold on;plot(s.m.p(eini(:),1),s.m.p(eini(:),2),'*')
hold on;plot(s.m.p(eout(:),1),s.m.p(eout(:),2),'*')
axis off equal


s.m = s.m.evalKeFeC('TL6');
%%
statorLoss = 32.6;
copperLoss = 16.9;
hin = 78;
hout = 5.5;
hend = 5.5;
Tairgap = 25;
Tambient = 25;
s = s.clearallbcs;
s = s.setrbc([eini;einc],hin*1e-3,Tairgap);
s = s.setrbc(eout,hout*1e-3,Tambient);
s = s.setExcitation('s1',statorLoss/24*(1000/55),'L');
%%
hold on
Niter = 4;
% myT = zeros(4,3,Niter);
Tcopper = zeros(1,Niter+1);
Tcopper(1) = 30;

for i = 1:Niter
    tmp = (copperLoss*(1+3.81e-3*(Tcopper(i)-25)))/24 *(1000/55);
    s = s.setExcitation('c1',tmp,'L');
    s = s.assignEdata;
    s = s.solve;
    myT = s.getMMMTemp;
    Tcopper(i+1) = myT(3,1);
end
% s.plotT
% s.m.plotwf
plot(Tcopper,'marker','s')
set(gca,'xtick',1:Niter+1)
xlabel('Iteration')
ylabel('Average Winding Temperature [C]')
title('2D Heat Transfer Simulation')
%%
clc
(copperLoss*(1+3.81e-3*(Tcopper(end)-25)))/24/55/s.m.mzs.c1.area
statorLoss/24/55/s.m.mzs.s1.area