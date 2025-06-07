%% initialization
clc
clear
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
p.Np = 188;
p.gammar = asin(p.R1*sin(p.betar/2)/p.R0);
p.gammas = asin((p.R1+p.g)*sin(p.betas/2)/p.R2);
% mesh size in rotor, stator and air gap
p.hr = 2;
p.hs =2;
p.hg = 1.5;
p.thetag = 0.8;
%% creation of rotor
p1 = [p.rsh,0];
p2 = [p.R1,0];
p3 = p.R1*[cos(p.betar/2),sin(p.betar/2)];
p4 = p.R0*[cos(p.gammar),sin(p.gammar)];
p5 = p.R0*[cos(pi/8),sin(pi/8)];
p6 = p.rsh*[cos(pi/8),sin(pi/8)];
g = g.newdlinewdkps(p2,p1,'l1',min(p.hg,2*p.R1*sin(p.thetag*pi/360)),'l2',p.hr);
g = g.newdarccppwdkps([0,0],p2,p3,'maxLength',min(p.hr,p.hg),'maxDegree',p.thetag);
g = g.newdlinewdkps(p3,p4,'l1',2*p.R1*sin(p.thetag*pi/360),'l2',2*p.R0*sin(3*pi/360));
g = g.newdarccppwdkps([0,0],p4,p5,'maxDegree',3);
g = g.newdlinewdkps(p5,p6,'l1',2*p.R0*sin(3*pi/360),'l2',p.hr);
g = g.newdarccppwdkps([0,0],p6,p1,'direction',-1,'maxLength',p.hr,'maxDegree',15);
g = g.filletla('L2','A2',2.5,'Nnodes',4);
g = g.newcb('r1','L1',-1,'A1',1,'L2',1,'A4',-1,'A2',1,'L3',1,'A3',1);
%% creation of stator
p1 = [p.R1+p.g,0];
p2 = [p.R3,0];
p3 = p.R3*[cos(pi/12),sin(pi/12)];
p4 = p.R2*[cos(pi/12),sin(pi/12)];
p5 = p.R2*[cos(p.gammas),sin(p.gammas)];
p6 = (p.R1+p.g)*[cos(p.betas/2),sin(p.betas/2)];
g = g.newdlinewdkps(p1,p2,'l1',min(p.hg,2*(p.R1+p.g)*sin(p.thetag*pi/360)),'l2',p.hs);
g = g.newdarccppwdkps([0,0],p2,p3,'maxLength',p.hs,'maxDegree',5);
g = g.newdlinewdkps(p3,p4,'maxLength',0.95*p.hs);
g = g.newdarccppwdkps([0,0],p4,p5,'direction',-1,'maxLength',min(p.hs,p.hg),'maxDegree',4);
g = g.newdlinewdkps(p6,p5,'l1',min(p.hg,2*(p.R1+p.g)*sin(p.thetag*pi/360)),...
    'l2',min(p.hg,2*p.R2*sin(4*pi/360)));
g = g.newdarccppwdkps([0,0],p6,p1,'direction',-1,'maxLength',min(p.hs,p.hg),'maxDegree',p.thetag);
g = g.filletla('L6','A6',2.5,'Nnodes',5);
g = g.newcb('s1','L4',1,'A5',1,'L5',1,'A6',1,'A8',-1,'L6',-1,'A7',1);
%% creation of coil
pw = (p.R1+p.g)*[cos(pi/12),sin(pi/12)];
g = g.newdlinewdkps(pw,p4,'l1',min(p.hg,2*(p.R1+p.g)*sin(p.thetag*pi/360)),...
    'l2',min(p.hg,2*p.R2*sin(4*pi/360)));
g = g.newdarccppwdkps([0,0],pw,p6,'direction',-1,'maxDegree',p.thetag);
g = g.newcb('c1','L6',1,'A8',1,'A6',-1,'L7',-1,'A9',1);
%% creation of rotor air pockets
g = g.cmirrorbd('L2',[cos(pi/8),sin(pi/8)]);
g = g.cmirrorbd('A2',[cos(pi/8),sin(pi/8)]);
g = g.cmirrorbd('A4',[cos(pi/8),sin(pi/8)]);
g = g.newdarccpp('kp3','kp18',[0,0],'maxLength',0.95*p.hg,'maxDegree',p.thetag);
g = g.newcb('rap1','L2',-1,'A12',1,'L8',1,'A11',-1 ...
    ,'A10',1,'A2',-1,'A4',1);
g = g.newdDM('r1','r1');
g = g.newdDM('s1','s1');
g = g.newdDM('c1','c1');
g = g.newdDM('rap1','rap1');
close all
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

mz = TMZPC(m.mzs.stator.t,m.mzs.stator.p);
mz = mz.moveNodes;
m.mzs.stator.p = mz.nodes;

mz = TMZPC(m.mzs.rotor.t,m.mzs.rotor.p);
mz = mz.moveNodes;
m.mzs.rotor.p = mz.nodes;

m = m.ggmesh;
szi = m.mzs.stator.zi;
rzi = m.mzs.rotor.zi;
Nts = sum(m.tzi(:,szi));
Ntr = sum(m.tzi(:,rzi));
Nt = m.Ngt;
%% solver setting
s = IHNLNRMSTL3(m);clear m;
s.scs.l = 1e-3;
s.scs.f = 1e6;
thetao = 0;
%% proccess

rpm = 1500;
f = fopen('C:\Users\AliJamalifard\Desktop\ia.tab','r');
i_phaseA = fscanf(f,'%f');
fclose(f);
i_phaseA = reshape(i_phaseA,2,[]);
i_phaseA(1,:) = rpm*(pi/30)*i_phaseA(1,:)*1e-3;
simulationAngle = pi/2;

ali = tic;
Nsim = 10;
Niter = 7;
tmp = Nsim;
for i = 1:Niter
    tmp = tmp + (tmp-1);
end
xi = linspace(0,simulationAngle,Nsim);
Ntheta = length(xi);

% Bxs = zeros(Nts,Ntheta);
% Bys = zeros(Nts,Ntheta);
% Bxr = zeros(Ntr,Ntheta);
% Byr = zeros(Ntr,Ntheta);

Bx = zeros(Nt,tmp);
By = zeros(Nt,tmp);

aLFlux = zeros(Ntheta,1);
bLFlux = zeros(Ntheta,1);
cLFlux = zeros(Ntheta,1);


for i = 1:Ntheta
    
    s.m = s.m.rotatemz('rotor',xi(i)-thetao);
    s.m = s.m.rotatemz('rap1',xi(i)-thetao);
    s.m = s.m.rotatemz('rap2',xi(i)-thetao);
    thetao = xi(i);
    
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
    
    
    s.m = s.m.evalKeFeC('TL3');
    s = s.clearallbcs;
    s = s.setdbc(k0,0);
    s = s.setopbc(km,ks);
    
    ia = interp1(i_phaseA(1,:),i_phaseA(2,:),rem(xi(i),pi/4));
    ib = interp1(i_phaseA(1,:),i_phaseA(2,:),rem(xi(i)+15*pi/180,pi/4));
    ic = interp1(i_phaseA(1,:),i_phaseA(2,:),rem(xi(i)+30*pi/180,pi/4));
    
    s = s.setExcitation('c1',ia*188,'C');
    s = s.setExcitation('c6',ia*188,'C');
    s = s.setExcitation('c2',ib*188,'C');
    s = s.setExcitation('c3',-ib*188,'C');
    s = s.setExcitation('c4',-ic*188,'C');
    s = s.setExcitation('c5',ic*188,'C');
    
    s = s.assignEdata;
    s = s.solve(1e-6,10);

    Bx(:,i) = s.B(1:Nt,1);
    By(:,i) = s.B(1:Nt,2);
    
    aLFlux(i) = s.evalLF('c1')+s.evalLF('c6');
    bLFlux(i) = s.evalLF('c2')-s.evalLF('c3');
    cLFlux(i) = s.evalLF('c5')-s.evalLF('c4');
    

    
end

Intervals = [1:Nsim-1;2:Nsim]';
Flags = ones(1,Nsim-1);

sLoss= zeros(1,Niter+1);
rLoss= zeros(1,Niter+1);
seLoss= zeros(1,Niter+1);
reLoss= zeros(1,Niter+1);

myNsim= zeros(1,Niter+1);
myIndex = Nsim;
[sLoss(1),rLoss(1)] = HysteresisLossGSEnew(Bx(:,1:Nsim),By(:,1:Nsim),s,xi);
[seLoss(1),reLoss(1)] = EdyyCurrentNew(Bx(:,1:Nsim),By(:,1:Nsim),s,xi);

myNsim(1) = Nsim;

for iter = 1:Niter
    
    n = size(Intervals,1);
    if ~any(Flags)
%         clc
        disp(iter)
        break
    end
        
    
    for i = 1:n
        
        if Flags(i)
            
            xtmp = mean(xi(Intervals(i,:)));
            ytmp = [mean(Bx(:,Intervals(i,:)),2),...
                mean(By(:,Intervals(i,:)),2)];
            
            s.m = s.m.rotatemz('rotor',xtmp-thetao);
            s.m = s.m.rotatemz('rap1',xtmp-thetao);
            s.m = s.m.rotatemz('rap2',xtmp-thetao);
            thetao = xtmp;
            
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
            
            
            s.m = s.m.evalKeFeC('TL3');
            s = s.clearallbcs;
            s = s.setdbc(k0,0);
            s = s.setopbc(km,ks);
            
            ia = interp1(i_phaseA(1,:),i_phaseA(2,:),rem(xtmp,pi/4));
            ib = interp1(i_phaseA(1,:),i_phaseA(2,:),rem(xtmp+15*pi/180,pi/4));
            ic = interp1(i_phaseA(1,:),i_phaseA(2,:),rem(xtmp+30*pi/180,pi/4));
            
            s = s.setExcitation('c1',ia*188,'C');
            s = s.setExcitation('c6',ia*188,'C');
            s = s.setExcitation('c2',ib*188,'C');
            s = s.setExcitation('c3',-ib*188,'C');
            s = s.setExcitation('c4',-ic*188,'C');
            s = s.setExcitation('c5',ic*188,'C');
            
            s = s.assignEdata;
            s = s.solve(1e-6,10);
    
            bmid = s.B(1:Nt,:);
            bmid =  sqrt(sum((ytmp-bmid).^2,2));
            
            mznames = fieldnames(s.m.mzs);
            Nu = ones(1,s.m.Ngt);
            for ii = 1:s.m.Nmzs-1
                    mzname = mznames{ii};
                    if ~s.m.mts.(s.m.mzs.(mzname).material).MagneticPermeability.isLinear
                        Nu(s.m.tzi(:,s.m.mzs.(mzname).zi)) =...
                            ppval(s.m.mts.(s.m.mzs.(mzname).material).vB,...
                            bmid(s.m.tzi(:,s.m.mzs.(mzname).zi)));
                    end
            end
                
errnorm = 0.5*1e-6*7.9577e+05*(Nu(1:Nt).*(s.m.gta(1:Nt)))* bmid.^2;

% errnorm =  s.m.gta(1:Nt)' .* sum(abs(ytmp-bmid),2) /max(s.m.gta(1:Nt));

            myIndex = myIndex+1;
%             mean(tmp2./tmp1)>0.01
% errnorm>1e-4
            if errnorm>1
                xi = [xi,xtmp];
                Bx(:,myIndex) = s.B(1:Nt,1);
                By(:,myIndex) = s.B(1:Nt,2);

                aLFlux(end+1) = s.evalLF('c1')+s.evalLF('c6');
                bLFlux(end+1) = s.evalLF('c2')-s.evalLF('c3');
                cLFlux(end+1) = s.evalLF('c5')-s.evalLF('c4');
                Intervals(end+1,:) = [length(xi),Intervals(i,2)];
                Intervals(i,2) = length(xi);
                Flags(end+1) = 1;
            else
                xi = [xi,xtmp];
                Bx(:,myIndex) = s.B(1:Nt,1);
                By(:,myIndex) = s.B(1:Nt,2);

                aLFlux(end+1) = s.evalLF('c1')+s.evalLF('c6');
                bLFlux(end+1) = s.evalLF('c2')-s.evalLF('c3');
                cLFlux(end+1) = s.evalLF('c5')-s.evalLF('c4');
                Intervals(end+1,:) = [length(xi),Intervals(i,2)];
                Intervals(i,2) = length(xi);
                Flags(end+1) = 0;
                Flags(i) = 0;
            end
        end
        
    end
    
    [~,index] = sort(xi);
%     ttmp = index(1:Nsim)
    [sLoss(1+iter),rLoss(1+iter)] = HysteresisLossGSEnew(Bx(:,index),By(:,index),s,sort(xi));
    [seLoss(1+iter),reLoss(1+iter)] = EdyyCurrentNew(Bx(:,index),By(:,index),s,sort(xi));
    
    myNsim(1+iter) = length(xi);
end

clc
toc(ali)

Bx = Bx(:,1:myIndex);
By = By(:,1:myIndex);
[xi,index] = sort(xi);
Bx = Bx(:,index);
By = By(:,index);
% Bxs = Bxs(:,index);
% Bys = Bys(:,index);
% Bxr = Bxr(:,index);
% Byr = Byr(:,index);
aLFlux = aLFlux(index);
bLFlux = bLFlux(index);
cLFlux = cLFlux(index);
xrotorAngle = xi;

aLFlux = aLFlux*2*188*55/1000;
bLFlux = bLFlux*2*188*55/1000;
cLFlux = cLFlux*2*188*55/1000;
% 
close all
% subplot(211)
% hold all
% plot(xrotorAngle*180/pi,ia)
% plot(xrotorAngle*180/pi,ib)
% plot(xrotorAngle*180/pi,ic)
% xlabel('Rotor Angle [Deg]')
% ylabel('Phase Current [A]')
% legend('Phase A','Phase B','Phase C')
% subplot(212)
hold all
plot(xrotorAngle*180/pi,aLFlux,'linewidth',1.5)
plot(xrotorAngle*180/pi,bLFlux,'linewidth',1.5)
plot(xrotorAngle*180/pi,cLFlux,'linewidth',1.5)
xlabel('Rotor Angle [Deg]')
ylabel('Phase Linkage Flux [wb]')
title('Adaptive Method')
legend('Phase A','Phase B','Phase C')
stem(xrotorAngle*180/pi,max([aLFlux';bLFlux';cLFlux']),'filled','color',[136,136,136]/255)
set(gca,'xlim',[0,45])
% 
% Torque = Torque*(55/0.3/4/pi/100);