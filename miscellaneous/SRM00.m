%% initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
g = G2DKERNEL;
% parameters
p = ID_SRM0(18,12,[cd,'\srm2_8_3.win']);
p.D = 0.65;
p.alphas = 0.5;
p.alphar = 0.42;
p.hrPg = 25;
% mesh size in rotor, stator and air gap
thetar1 = 8*pi/180;
thetar2 = 5*pi/180;
thetag = 1.5*pi/180;
thetas1 = 3*pi/180;
thetas2 = 3*pi/180;

msr1 = 2*p.rsh*sin(thetar1/2);
msr2 = 2*p.Rrm*sin(thetar2/2);
msg = 1.5*(p.Rro+p.g/2)*sin(thetag/2);
mss1 = 2*p.Rsm*sin(thetas1/2);
mss2 = 2*p.Rso*sin(thetas2/2);
%% creation of rotor
p1 = [p.rsh,0];
p2 = [p.Rro,0];
p3 = p.Rro*[cos(p.betar/2),sin(p.betar/2)];
p4 = p.Rrm*[cos(p.gammar),sin(p.gammar)];
p5 = p.Rrm*[cos(p.Tr/2),sin(p.Tr/2)];
p6 = p.rsh*[cos(p.Tr/2),sin(p.Tr/2)];
g = g.newdswdkps(p1,p2,'l1',msr1,'l2',msg);
g = g.newdarccppwdkps([0,0],p2,p3,'maxLength',msg);
g = g.newdswdkps(p3,p4,'l1',msg,'l2',msr2);
g = g.newdarccppwdkps([0,0],p4,p5,'maxLength',msr2);
g = g.newdswdkps(p5,p6,'l1',msr2,'l2',msr1);
g = g.newdarccppwdkps([0,0],p6,p1,'direction',-1,'maxLength',msr1);
g = g.newloop('r1','S1',1,'A1',1,'S2',1,'A2',1,'S3',1,'A3',1);
%% creation of stator
p1 = [p.Rro+p.g,0];
p2 = [p.Rso,0];
p3 = p.Rso*[cos(p.Ts/2),sin(p.Ts/2)];
p4 = p.Rsm*[cos(p.Ts/2),sin(p.Ts/2)];
p5 = p.Rsm*[cos(p.gammas),sin(p.gammas)];
p6 = (p.Rro+p.g)*[cos(p.betas/2),sin(p.betas/2)];
g = g.newdswdkps(p1,p2,'l1',msg,'l2',mss2);
g = g.newdarccppwdkps([0,0],p2,p3,'maxLength',mss2);
g = g.newdswdkps(p3,p4,'l1',mss2,'l2',mss1);
g = g.newdarccppwdkps([0,0],p4,p5,'direction',-1,'maxLength',mss1);
g = g.newdswdkps(p6,p5,'l1',msg,'l2',mss1);
g = g.newdarccppwdkps([0,0],p6,p1,'direction',-1,'maxLength',msg);
g = g.newloop('s1','S4',1,'A4',1,'S5',1,'A5',1,'S6',-1,'A6',1);
%% creation of coil
pw = (p.Rro+p.g)*[cos(p.Ts/2),sin(p.Ts/2)];
g = g.newdswdkps(pw,p4,'l1',msg,'l2',mss1);
g = g.newdarccppwdkps([0,0],pw,p6,'direction',-1,'maxLength',msg);
g = g.newloop('c1','S6',1,'A5',-1,'S7',-1,'A7',1);
%% creation of rotor air pockets
g = g.cmirrored('S2',[cos(p.Tr/2),sin(p.Tr/2)]);
g = g.cmirrored('A2',[cos(p.Tr/2),sin(p.Tr/2)]);
g = g.newdarccpp('kp3','kp14',[0,0],'maxLength',msg);
g = g.newloop('rap1','S2',-1,'A9',1,'S8',1,'A8',1,'A2',-1);
%% creation of domains
g = g.newface('r1','r1');
g = g.newface('s1','s1');
g = g.newface('c11','c1');
g = g.newface('rap1','rap1');
%% creation of mesh
m = TMDBC();
m = m.addmz('s1',g.getdmmz('s1'));
m = m.addmz('r1',g.getdmmz('r1'));
m = m.addmz('c11',g.getdmmz('c11'));
m = m.addmz('rap1',g.getdmmz('rap1'));
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
m = m.setmzColor('c21',[0,181,250]/256);
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


m = m.ggmesh;
m.showmzs;
kr = m.getnIndexOnCircle([0,0],p.Rro);
ks = m.getnIndexOnCircle([0,0],p.Rro+p.g);
rps = m.nodes(kr,:);
sps = m.nodes(ks,:);

s = IHNLNRMSTL3(m);clear m
s.units.length = 'mm';
s.units.currentDensity = 'A/mm^2';
s.units.magneticVectorPotential = 'A/m';

%% proccess
phaseCurrent = linspace(0,10,5);
rotorAngle = linspace(0,p.Tr/2,5);
Ni = length(phaseCurrent);
Ntheta = length(rotorAngle);
linkageFlux = zeros(Ntheta,Ni);
theta_old = 0;
% loop for calculation of linkage flux 
for i = 1:Ntheta
    % rotation of regions
    s.m = s.m.rotatemz('rotor',rotorAngle(i)-theta_old);
        for j = 1:p.nr/p.Nsplit
            s.m = s.m.rotatemz(['rap',num2str(j)],rotorAngle(i)-theta_old);
        end
    if p.Nsplit>1
        s = s.addmz('AG',getmz(rotaterps(MC_AAG(sps,rps),rotorAngle(i))));
    else
        s = s.addmz('AG',getmz(rotaterps(MC_AG(sps,rps),rotorAngle(i))));
    end
    theta_old = rotorAngle(i);
    s.m = s.m.ggmesh;
    % getting index fo boundary conditions
    k0 = [s.m.getnIndexOnCircle([0,0],p.rsh);...
        s.m.getnIndexOnCircle([0,0],p.Rso)];
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
    s = s.removemz('AG');
end
linkageFlux = linkageFlux*p.nser*p.Lst*p.Np/p.np/1000;
close all
hold all
tmpC = [phaseCurrent,fliplr(phaseCurrent)];
tmpL = [linkageFlux(1,:),fliplr(linkageFlux(end,:))];
for i=1:Ntheta
    for j = 1:Ni
        plot(phaseCurrent,linkageFlux(i,:),...
            'color','b','marker','s','markersize',8,'linestyle','--',...
            'linewidth',1.5);
    end
end
title(['Wcon = ',num2str(polyarea(tmpC,tmpL)),' Wcon,des = ',num2str(p.Wcon)]);
xlabel('Phase Current [A]');
ylabel('Linkage Flux [wb]');
