%% initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
g = GDBC2D;
% parameters
p = ID_SRM0(18,12,[cd,'\srm18_12_3.win']);
p.D = 0.65;
p.alphas = 0.4;
p.alphar = 0.44;
p.hrPg = 25;
% mesh size in rotor, stator and air gap
thetar1 = 6*pi/180;
msr1 = 2*p.rsh*sin(thetar1/2);
thetar2 = 3.5*pi/180;
msr2 = 2*p.Rrm*sin(thetar2/2);
thetag = 1.5*pi/180;
msg = 2*(p.Rro+p.g/2)*sin(thetag/2);
thetas1 = 2.5*pi/180;
mss1 = 2*p.Rsm*sin(thetas1/2);
thetas2 = 4*pi/180;
mss2 = 2*p.Rso*sin(thetas2/2);
%% creation of rotor
p1 = [p.rsh,0];
p2 = [p.Rro,0];
p3 = p.Rro*[cos(p.betar/2),sin(p.betar/2)];
p4 = p.Rrm*[cos(p.gammar),sin(p.gammar)];
p5 = p.Rrm*[cos(p.Tr/2),sin(p.Tr/2)];
p6 = p.rsh*[cos(p.Tr/2),sin(p.Tr/2)];
g = g.newdlinewdkps(p1,p2,'l1',msr1,'l2',msg);
g = g.newdarccppwdkps([0,0],p2,p3,'maxLength',msg);
g = g.newdlinewdkps(p3,p4,'l1',msg,'l2',msr2);
g = g.newdarccppwdkps([0,0],p4,p5,'maxLength',msr2);
g = g.newdlinewdkps(p5,p6,'l1',msr2,'l2',msr1);
g = g.newdarccppwdkps([0,0],p6,p1,'direction',-1,'maxLength',msr1);
g = g.newcb('r1','L1',1,'A1',1,'L2',1,'A2',1,'L3',1,'A3',1);
%% creation of stator
p1 = [p.Rro+p.g,0];
p2 = [p.Rso,0];
p3 = p.Rso*[cos(p.Ts/2),sin(p.Ts/2)];
p4 = p.Rsm*[cos(p.Ts/2),sin(p.Ts/2)];
p5 = p.Rsm*[cos(p.gammas),sin(p.gammas)];
p6 = (p.Rro+p.g)*[cos(p.betas/2),sin(p.betas/2)];
g = g.newdlinewdkps(p1,p2,'l1',msg,'l2',mss2);
g = g.newdarccppwdkps([0,0],p2,p3,'maxLength',mss2);
g = g.newdlinewdkps(p3,p4,'l1',mss2,'l2',mss1);
g = g.newdarccppwdkps([0,0],p4,p5,'direction',-1,'maxLength',mss1);
g = g.newdlinewdkps(p6,p5,'l1',msg,'l2',mss1);
g = g.newdarccppwdkps([0,0],p6,p1,'direction',-1,'maxLength',msg);
g = g.newcb('s1','L4',1,'A4',1,'L5',1,'A5',1,'L6',-1,'A6',1);
%% creation of coil
pw = (p.Rro+p.g)*[cos(p.Ts/2),sin(p.Ts/2)];
g = g.newdlinewdkps(pw,p4,'l1',msg,'l2',mss1);
g = g.newdarccppwdkps([0,0],pw,p6,'direction',-1,'maxLength',msg);
g = g.newcb('c1','L6',1,'A5',-1,'L7',-1,'A7',1);
%% creation of rotor air pockets
g = g.cmirrorbd('L2',[cos(p.Tr/2),sin(p.Tr/2)]);
g = g.cmirrorbd('A2',[cos(p.Tr/2),sin(p.Tr/2)]);
g = g.newdarccpp('kp3','kp14',[0,0],'maxLength',msg);
g = g.newcb('rap1','L2',-1,'A9',1,'L8',1,'A8',1,'A2',-1);
%% creation of domains
g = g.newdDM('r1','r1');
g = g.newdDM('s1','s1');
g = g.newdDM('c11','c1');
g = g.newdDM('rap1','rap1');
close all
%% creation of mesh
m = MDBCT(g);clear g;
m = m.setmzColor('r1','r');
m = m.setmzColor('s1','r');
m = m.setmzColor('c11','m');
m = m.setmzColor('rap1','w');
m = m.setMaterial('r1','m19_24ga');
m = m.setMaterial('s1','m19_24ga');
m = m.setMaterial('c11','copper');
m = m.addMaterial(materialdir,'air');
m = m.addMaterial(materialdir,'m19_24ga');
m = m.addMaterial(materialdir,'copper');
m = m.cmirrormz('r2','r1',[cos(p.Tr/2),sin(p.Tr/2)]);
m = m.cmirrormz('s2','s1',[1,0]);
m = m.cmirrormz('c21','c11',[1,0]);
for i = 1:2:2*(p.nr/p.Nsplit-1)
    m = m.crotatemz(['r',num2str(i+2)],['r',num2str(i)],p.Tr);
    m = m.crotatemz(['r',num2str(i+3)],['r',num2str(i+1)],p.Tr);
end
for i = 1:2:2*(p.ns/p.Nsplit-1)
    m = m.crotatemz(['s',num2str(i+2)],['s',num2str(i)],p.Ts);
    m = m.crotatemz(['s',num2str(i+3)],['s',num2str(i+1)],p.Ts);
end
for i = 1:p.ns/p.Nsplit-1
    m = m.crotatemz(['c1',num2str(i+1)],['c1',num2str(i)],p.Ts);
    m = m.crotatemz(['c2',num2str(i+1)],['c2',num2str(i)],p.Ts);
end
for i = 1:p.nr/p.Nsplit-1
    m = m.crotatemz(['rap',num2str(i+1)],['rap',num2str(i)],p.Tr);
end
% joining meshs
temp = getlist('r',1:2*p.nr/p.Nsplit);
m = m.joinmzs('rotor',temp{:});
temp = getlist('s',1:2*p.ns/p.Nsplit);
m = m.joinmzs('stator',temp{:});
m = m.rotatemz('rotor',-p.Tr/2);
for i = 1:p.nr/p.Nsplit
    m = m.rotatemz(['rap',num2str(i)],-p.Tr/2);
end
%% solver setting
ali = tic;
s = IHNLNRMSTL3old(m);clear m;
% setting units
s.scs.l = 1e-3;
s.scs.f = 1e6;
%% proccess
phaseCurrent = 0:10;
rotorAngle = linspace(0,p.Tr/2,10);
Ni = length(phaseCurrent);
Ntheta = length(rotorAngle);
linkageFlux = zeros(Ntheta,Ni);
% loop for calculation of linkage flux 
for i = 1:Ntheta
    s.m = s.m.ggmesh;
    kin = s.m.getIndexOnCircle([0,0],p.Rro);
    kout = s.m.getIndexOnCircle([0,0],p.Rro+p.g);
    if p.Nsplit>1
        s.m = s.m.makeAAG(kin,kout);
    else
        s.m = s.m.makeAG(kin,kout);
    end
    % getting index fo boundary conditions
    k0 = [s.m.getIndexOnCircle([0,0],p.rsh);...
        s.m.getIndexOnCircle([0,0],p.Rso)];
    s = s.clearallbcs;
    s = s.setdbc(k0,0);
    if p.Nsplit>1
        k = s.m.getfb;
        k = setdiff(unique(k(:)),k0);
        [km,ks] = s.m.splitPeriodic(k,2*pi/p.Nsplit);
        s = s.setopbc(km,ks);
    end
    % calling solver
    s.m = s.m.evalKeFeC('TL3');
    for j = 1:Ni
        % set excitation
        for k = 1:size(p.win.phase1,1)
            if p.win.phase1(k,1) < p.ns/p.Nsplit
                if p.win.phase1(k,2) > 0
                    s = s.setExcitation(['c1',num2str(p.win.phase1(k,1))],...
                        phaseCurrent(j)*p.Np/p.np,'C');
                    s = s.setExcitation(['c2',num2str(p.win.phase1(k,1))],...
                        -phaseCurrent(j)*p.Np/p.np,'C');
                else
                    s = s.setExcitation(['c1',num2str(p.win.phase1(k,1))],...
                        -phaseCurrent(j)*p.Np/p.np,'C');
                    s = s.setExcitation(['c2',num2str(p.win.phase1(k,1))],...
                        phaseCurrent(j)*p.Np/p.np,'C');
                end
            end
        end
        s = s.assignEdata;
        s = s.solve(1e-8,30);
        % eval linkage flux
        for k = 1:size(p.win.phase1,1)
            if p.win.phase1(k,1) < p.ns/p.Nsplit
                if p.win.phase1(k,2) > 0
                    linkageFlux(i,j) = linkageFlux(i,j) ...
                        + s.evalLF(['c1',num2str(p.win.phase1(k,1))]) ...
                        - s.evalLF(['c2',num2str(p.win.phase1(k,1))]);
                else
                    linkageFlux(i,j) = linkageFlux(i,j)  ...
                        - s.evalLF(['c1',num2str(p.win.phase1(k,1))]) ...
                        + s.evalLF(['c2',num2str(p.win.phase1(k,1))]);
                end
            end
        end
    end
    if i<Ntheta
        s.m = s.m.rotatemz('rotor',rotorAngle(2));
        for j = 1:p.nr/p.Nsplit
            s.m = s.m.rotatemz(['rap',num2str(j)],rotorAngle(2));
        end
    end 
end
clc
toc(ali)
linkageFlux = linkageFlux*p.nser*p.Lst*p.Np/p.np/1000;
close all
hold all
tmpC = [phaseCurrent,fliplr(phaseCurrent)];
tmpL = [linkageFlux(1,:),fliplr(linkageFlux(end,:))];
for i=1:Ntheta
    for j = 1:Ni
        plot(phaseCurrent*2,linkageFlux(i,:),...
            'color','b','marker','s','markersize',8,'linestyle','--',...
            'linewidth',1.5);
    end
end
title(['Wcon = ',num2str(polyarea(tmpC,tmpL)),' Wcon,des = ',num2str(p.Wcon)]);
xlabel('Phase Current [A]');
ylabel('Linkage Flux [wb]');