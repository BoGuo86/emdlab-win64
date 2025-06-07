%% initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
g = GDBC2D;
%% permitted flux densities [Tesla]
p.Bg = 0.87;
p.Bt = 1.7;
p.Bsy = 1.4;

p.Bry = 1.4;

p.Kst = 0.97;
p.Kf = 0.45;

%% dimensions [mm]
% rated speed [rpm]
p.rpm = 500;
% rated mechanical angular speed
p.wm = p.rpm*2*pi/60;
% inner stator bore diameter
p.Rso = 50;
% 
p.D = 0.6;
p.Rro = p.D*p.Rso;
% stack length
p.Lst = 100;
% number of stator slot
p.Ns = 24;
% number of poles
p.Nm = 8;
% magnet height
p.hm = 4;
% air gap length
p.g = 1;

p.wtb = (2*pi*p.Rro*p.Bg)/(p.Ns*p.Kst*p.Bt);
p.wtb = round(p.wtb,1);
p.wsy = (pi*p.Rro*p.Bg)/(p.Nm*p.Kst*p.Bsy);
p.wsy = round(p.wsy,1);
p.wry = (pi*p.Rro*p.Bg)/(p.Nm*p.Kst*p.Bry);
p.wry = round(p.wry,1);

p.rsh = p.Rro - p.hm - p.wry;

% % getting balanced 3phase winding layout
% [d,phaseA,phaseB,phaseC] = Winding(p.Nm,p.Ns);
% angular pitch of stator slot
p.Ts = 2*pi/p.Ns;
% slot opening angles
p.alphas = 0.3;
p.Tso1 = p.alphas*p.Ts;

p.Tso2 = p.Ts - 2*asin((p.wtb/2)/(p.Rso-p.wsy));
% depth of slot
p.hss1 = 1.5;
p.hss2 = 0.5;
p.hss3 = p.Rso - p.wsy - p.Rro - p.g - p.hss1 - p.hss1*cos((p.Ts-p.Tso1)/2);

% mesh sizes
p.smsize = 1;
p.gmsize = 1;
p.rmsize = 1;
p.thetag = 1;
%% creation of stator lamination
p1 = [p.Rro+p.g,0];
p2 = [p.Rso,0];
p3 = pmove(p2,'theta',p.Ts/2);
p4 = pmove(p3,'r',-p.wsy);
p5 = pmove(p4,'theta',-p.Tso2/2);
p6 = pmove(p5,'x',-p.hss3);
p8 = pmove(p1,'theta',(p.Ts-p.Tso1)/2);
p7 = pmove(p8,'r',p.hss1);
[g,L1] = g.newdlinewdkps(p1,p2,p.smsize);
[g,A1] = g.newdarccppwdkps([0,0],p2,p3,1,p.smsize,5);
[g,L2] = g.newdlinewdkps(p3,p4,p.smsize);
[g,A2] = g.newdarccppwdkps([0,0],p4,p5,-1,min(p.smsize,p.gmsize),5);
[g,L3] = g.newdlinewdkps(p5,p6,min(p.smsize,p.gmsize));
[g,L4] = g.newdlinewdkps(p6,p7,min(p.smsize,p.gmsize));
[g,L5] = g.newdlinewdkps(p7,p8,min(p.smsize,p.gmsize));
[g,A3] = g.newdarccppwdkps([0,0],p8,p1,-1,min(p.smsize,p.gmsize),p.thetag);
g = g.newcb('s1',L1,1,A1,1,L2,1,A2,1,L3,1,L4,1,L5,1,A3,1);
%% cration of coil
p9 = pmove(p1,'theta',p.Ts/2);
p9 = pmove(p9,'r',p.hss1);
g = g.newdlinewdkps(p7,p9,p.gmsize);
g = g.newdlinewdkps(p4,p9,p.gmsize);
g = g.newcb('c11','L3',-1,'A2',-1,'L7',1,'L6',-1,'L4',-1);
% creation of slot air
g = g.cmirrorbd('L5',[cos(p.Ts/2),sin(p.Ts/2)]);
g = g.cmirrorbd('L6',[cos(p.Ts/2),sin(p.Ts/2)]);
g = g.newdarccpp('kp8','kp11',[0,0],1,p.gmsize,p.thetag);
g = g.newcbwd('sap1',p.gmsize,'L5',-1,'L6',1,'L9',-1,'L8',1,'A4',-1);
%% magnet
p1 = [p.rsh+p.wry,0];
p2 = p1+[p.hm,0];
p3 = pmove(p2,'theta',2*pi/p.Nm/3);
p4 = pmove(p3,'r',-p.hm);
p5 = pmove(p4,'theta',1*pi/p.Nm/3);
p7 = [p.rsh,0];
p6 = pmove(p7,'theta',pi/p.Nm);
g = g.newdlinewdkps(p1,p2,p.gmsize);
g = g.newdarccppwdkps([0,0],p2,p3,1,p.gmsize,p.thetag);
g = g.newdlinewdkps(p3,p4,min(p.rmsize,p.gmsize));
g = g.newdarccppwdkps([0,0],p4,p1,-1,min(p.rmsize,p.gmsize),5);
g = g.newcb('m1','L10',1,'A5',1,'L11',1,'A6',1);
%% rotor lamination
g = g.newdarccppwdkps([0,0],p4,p5,1,min(p.rmsize,p.gmsize),5);
g = g.newdlinewdkps(p5,p6,p.rmsize);
g = g.newdarccppwdkps([0,0],p6,p7,-1,p.rmsize,5);
g = g.newdlinewdkps(p7,p1,p.rmsize);
g = g.newcb('r1','L13',1,'A6',-1,'A7',1,'L12',1,'A8',1);
%% rotor air pocket
g = g.cmirrorbd('L11',[cos(pi/p.Nm),sin(pi/p.Nm)]);
g = g.cmirrorbd('A7',[cos(pi/p.Nm),sin(pi/p.Nm)]);
g = g.newdarccpp('kp14','kp19',[0,0],1,p.gmsize,p.thetag);
g = g.newcb('rap1','L11',-1,'A10',1,'L14',1,'A9',1,'A7',-1);
%% creation of domains
g = g.newdDM('s1',p.smsize,'s1');
g = g.newdDM('c11',p.gmsize,'c11');
g = g.newdDM('m1',p.gmsize,'m1');
g = g.newdDM('r1',p.rmsize,'r1');
g = g.newdDM('rap1',p.gmsize,'rap1');
close all
%% creation of mesh
close all
m = MDBCT(g);clear g;
m = m.addMaterial(materialdir,'air');
m = m.addMaterial(materialdir,'m19_24ga');
m = m.addMaterial(materialdir,'copper');
m = m.setmzColor('s1',[0 191 255]/255);
m = m.setmzColor('r1',[0 191 255]/255);
m = m.setMaterial('s1','m19_24ga');
m = m.cmirrormz('s2','s1',[1,0]);
m = m.setMaterial('c11','copper');
m = m.setMaterial('r1','m19_24ga');
m = m.setmzColor('c11',[255 140 0]/255);
m = m.cmirrormz('c21','c11',[cos(p.Ts/2),sin(p.Ts/2)]);
m = m.setmzColor('c21','r');
m = m.setmzColor('sap1','w');
m = m.setmzColor('rap1','w');
for i = 1:2:2*(p.Ns-1)
    m = m.crotatemz(['s',num2str(i+2)],['s',num2str(i)],2*pi/p.Ns);
    m = m.crotatemz(['s',num2str(i+3)],['s',num2str(i+1)],2*pi/p.Ns);
end
m = m.cmirrormz('r2','r1',[1,0]);
m = m.cmirrormz('m2','m1',[1,0]);
m = m.joinmzs('magnet1','m1','m2');
m = m.setmzColor('magnet1','m');
for i = 1:(p.Ns-1)
    m = m.crotatemz(['c1',num2str(i+1)],['c1',num2str(i)],2*pi/p.Ns);
    m = m.crotatemz(['c2',num2str(i+1)],['c2',num2str(i)],2*pi/p.Ns);
    m = m.crotatemz(['sap',num2str(i+1)],['sap',num2str(i)],2*pi/p.Ns);
end
for i = 1:(p.Nm-1)
    m = m.crotatemz(['magnet',num2str(i+1)],['magnet',num2str(i)],2*pi/p.Nm);
    m = m.crotatemz(['rap',num2str(i+1)],['rap',num2str(i)],2*pi/p.Nm);
end
for i = 1:2:2*(p.Nm-1)
    m = m.crotatemz(['r',num2str(i+2)],['r',num2str(i)],2*pi/p.Nm);
    m = m.crotatemz(['r',num2str(i+3)],['r',num2str(i+1)],2*pi/p.Nm);
end
temp = getlist('s',1:2*p.Ns);
m = m.joinmzs('stator',temp{:});
temp = getlist('r',1:2*p.Nm);
m = m.joinmzs('rotor',temp{:});
m = m.ggmesh;
m.showmzs
szi = m.mzs.stator.zi;
rzi = m.mzs.rotor.zi;
Nts = sum(m.tzi(:,szi));
Ntr = sum(m.tzi(:,rzi));
% %% creation of solver
% s = IHNLNRMSTL3(m);clear m;
% s.scs.l = 1e-3;
% s.scs.f = 1e6;
% phaseA = phaseA +1;
% phaseB = phaseB +1;
% phaseC = phaseC +1;
% %%
% ia = 0;
% ib = 0;
% ic = 10;
% 
% 
% for j = 1:length(phaseA)
%     if d(j)>0
%         s = s.setExcitation(['c1',num2str(phaseA(j,1))],ia);
%         s = s.setExcitation(['c2',num2str(phaseA(j,2))],-ia);
%         s = s.setExcitation(['c1',num2str(phaseB(j,1))],ib);
%         s = s.setExcitation(['c2',num2str(phaseB(j,2))],-ib);
%         s = s.setExcitation(['c1',num2str(phaseC(j,1))],ic);
%         s = s.setExcitation(['c2',num2str(phaseC(j,2))],-ic);
%     else
%         s = s.setExcitation(['c1',num2str(phaseA(j,1))],-ia);
%         s = s.setExcitation(['c2',num2str(phaseA(j,2))],ia);
%         s = s.setExcitation(['c1',num2str(phaseB(j,1))],-ib);
%         s = s.setExcitation(['c2',num2str(phaseB(j,2))],ib);
%         s = s.setExcitation(['c1',num2str(phaseC(j,1))],-ic);
%         s = s.setExcitation(['c2',num2str(phaseC(j,2))],ic);
%     end
% end
% 
% s.m = s.m.ggmesh;
% k1 = s.m.getIndexOnCircle([0,0],p.rsh+p.bry+p.hm);
% k2 = s.m.getIndexOnCircle([0,0],p.rsh+p.bry+p.hm+p.g);
% s.m = s.m.makeAG(k1,k2);
% k0 = s.m.getfb;
% k0 = unique(k0(:));
% 
% s = s.clearallbcs;
% s = s.setdbc(k0,0);
% s.m = s.m.evalKeFeC('TL3');
% s = s.assignEdata;
% s = s.solve(1e-6,10);
% 
% aLFlux = 0;
% bLFlux = 0;
% cLFlux = 0;
% 
% for j = 1:length(d)
%     if d(j)>0
%         aLFlux = aLFlux + s.evalLF(['c1',num2str(phaseA(j,1))]);
%         aLFlux = aLFlux - s.evalLF(['c2',num2str(phaseA(j,2))]);
%         bLFlux = bLFlux + s.evalLF(['c1',num2str(phaseB(j,1))]);
%         bLFlux = bLFlux - s.evalLF(['c2',num2str(phaseB(j,2))]);
%         cLFlux = cLFlux + s.evalLF(['c1',num2str(phaseC(j,1))]);
%         cLFlux = cLFlux - s.evalLF(['c2',num2str(phaseC(j,2))]);
%     else
%         aLFlux = aLFlux - s.evalLF(['c1',num2str(phaseA(j,1))]);
%         aLFlux = aLFlux + s.evalLF(['c2',num2str(phaseA(j,2))]);
%         bLFlux = bLFlux - s.evalLF(['c1',num2str(phaseB(j,1))]);
%         bLFlux = bLFlux + s.evalLF(['c2',num2str(phaseB(j,2))]);
%         cLFlux = cLFlux - s.evalLF(['c1',num2str(phaseC(j,1))]);
%         cLFlux = cLFlux + s.evalLF(['c2',num2str(phaseC(j,2))]);
%     end
% end
% 
% 
% aLFlux = aLFlux*1000*p.Ncoil/p.Lstk
% bLFlux = bLFlux*1000*p.Ncoil/p.Lstk
% cLFlux = cLFlux*1000*p.Ncoil/p.Lstk
% 
% 
% s.plotBmag
% 
