%% initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
g = GDBC2D;
% parameters
p.rsh = 15.000000;
p.R0 = 30.400000;
p.R1 = 40.900000;
p.R2 = 65.800000;
p.R3 = 74.900000;
p.betar = 19 * pi/180;
p.betas = 16 * pi/180;
p.g = 0.3;
p.dw = 2;
p.Rw = p.R1+p.g+p.dw;
p.gammar = asin(p.R1*sin(p.betar/2)/p.R0);
p.gammas = asin((p.R1+p.g)*sin(p.betas/2)/p.R2);
p.gammaw = asin((p.R1+p.g)*sin(p.betas/2)/p.Rw);
% mesh size
p.hr = 3;
p.hs = 3;
p.hg = 1.5;
p.thetag = 1;
p.alpha = 3;
%% creation of rotor
p1 = pmove([p.rsh,0],'theta',-pi/8);
p2 = pmove([p.R0,0],'theta',-pi/8);
p3 = p.R0*[cos(-p.gammar),sin(-p.gammar)];
p4 = p.R1*[cos(-p.betar/2),sin(-p.betar/2)];
g = g.newdlinewdkps(p1,p2,p.hr);
g = g.newdarccppwdkps([0,0],p2,p3,1,min(p.hr,p.hg));
g = g.newdlinewdkps(p3,p4,p.hr);
g = g.cmirrorbd('L1',[1,0]);
g = g.cmirrorbd('A1',[1,0]);
g = g.newdarccpp('kp1','kp5',[0,0],1,p.hr,10);
pc = p.alpha*[cos(-pi/8),sin(-pi/8)];
p5 = g.findRayCircleIntersect(g.kps.kp7,[1 0],pc,norm(pc-g.kps.kp4));
g = g.newdlinewdkps(g.kps.kp7,p5,min(p.hr,p.hg));
g = g.newdarccpp('kp4','kp8',pc,1,min(p.hr,p.hg),p.thetag);
g = g.newcb('r1','L1',1,'A1',1,'L2',1,'A4',1,'L4',-1,'A2',-1,'L3',-1,'A3',-1);
g = g.newdDM('r1',p.hr,'r1');
%% creation ot rotor air gap
p6 = pmove([p.R1+p.g/3,0],'theta',-pi/8);
p7 = pmove([p.R1+p.g/3,0],'theta',pi/8);
g = g.newdlinewdkps(p2,p6,p.hg);
g = g.newdarccppwdkps([0,0],p6,p7,1,p.hg,p.thetag);
g = g.newdline('kp6','kp10',p.hg);
g = g.newcbwd('ra1',p.hg,'L5',1,'A5',1,'L6',-1,'A2',1,'L4',1,'A4',-1,'L2',-1,'A1',-1);
%% creation of stator
p1 = [p.R1+p.g,0];
p2 = [p.R3,0];
p3 = p.R3*[cos(pi/12),sin(pi/12)];
p4 = p.R2*[cos(pi/12),sin(pi/12)];
p5 = p.R2*[cos(p.gammas),sin(p.gammas)];
p6 = p.Rw*[cos(p.gammaw),sin(p.gammaw)];
p7 = (p.R1+p.g)*[cos(p.betas/2),sin(p.betas/2)];
g = g.newdlinewdkps(p1,p2,p.hs);
g = g.newdarccppwdkps([0,0],p2,p3,1,p.hs,5);
g = g.newdlinewdkps(p3,p4,p.hs);
g = g.newdarccppwdkps([0,0],p4,p5,-1,min(p.hs,p.hg),5);
g = g.newdlinewdkps(p5,p6,min(p.hs,p.hg));
g = g.newdlinewdkps(p6,p7,min(p.hs,p.hg));
g = g.newdarccppwdkps([0,0],p7,p1,-1,min(p.hs,p.hg),p.thetag);
g = g.newcbwd('s1',p.hs,'L7',1,'A6',1,'L8',1,'A7',1,'L9',1,'L10',1,'A8',1);
%% creation of coil
pw = p.Rw*[cos(pi/12),sin(pi/12)];
g = g.newdlinewdkps(p4,pw,p.hg);
g = g.newdarccppwdkps([0,0],pw,p6,-1,p.hg,5);
g =g.newcbwd('c1',p.hg,'L9',-1,'A7',-1,'L11',1,'A9',1);
%% creation of stator air gap
g = g.cmirrorbd('A9',[cos(pi/12),sin(pi/12)]);
g = g.cmirrorbd('L10',[cos(pi/12),sin(pi/12)]);
g = g.newdarccpp('kp17','kp20',[0,0],1,p.hg,p.thetag);
g =g.newcbwd('sa1',p.hg,'L10',-1,'A9',-1,'A10',1,'L12',1,'A11',-1);

%% meshing
m = MDBCT(g);clear g;
m = m.setMaterial('s1','m19_24ga');
m = m.setMaterial('r1','m19_24ga');
m = m.setMaterial('c1','copper');
m = m.addMaterial(materialdir,'air');
m = m.addMaterial(materialdir,'m19_24ga');
m = m.addMaterial(materialdir,'copper');
m = m.cmirrormz('s2','s1',[cos(pi/12),sin(pi/12)]);
m = m.cmirrormz('c2','c1',[1 0]);
for i = 1:2:22
    m = m.crotatemz(['s',num2str(i+2)],['s',num2str(i)],pi/6);
    m = m.crotatemz(['s',num2str(i+3)],['s',num2str(i+1)],pi/6);
    m = m.crotatemz(['c',num2str(i+2)],['c',num2str(i)],pi/6);
    m = m.crotatemz(['c',num2str(i+3)],['c',num2str(i+1)],pi/6);
end

for i = 2:12
    m = m.crotatemz(['sa',num2str(i)],['sa',num2str(i-1)],pi/6);
end

for i = 2:8
    m = m.crotatemz(['r',num2str(i)],['r',num2str(i-1)],pi/4);
    m = m.crotatemz(['ra',num2str(i)],['ra',num2str(i-1)],pi/4);
end


temp = getlist('s',1:24);
m = m.joinmzs('stator',temp{:});
temp = getlist('r',1:8);
m = m.joinmzs('rotor',temp{:});
temp = getlist('ra',1:8);
m = m.joinmzs('rair',temp{:});

% m.showmzs
% m.showmz('rotor')
% m.showmz('stator')
%% proccess
s = IHNLNRMSTL3(m);clear m;
s.scs.l = 1e-3;
s.scs.f = 1e6;
Np = 188;
%% solvings
cc = 5;
tt = linspace(0,45,50);
tt = [0,diff(tt)];
F(length(tt)) = struct('cdata',[],'colormap',[]);
w = zeros(length(tt),length(cc));
aa = zeros(length(tt),length(cc));
T = zeros(length(tt),length(cc));

for i = 1:length(tt)
    temp = tt(i)*pi/180;
    s.m = s.m.rotatemz('rotor',temp);
    s.m = s.m.rotatemz('rair',temp);
    s.m = s.m.ggmesh;
    k0 = [s.m.getIndexOnCircle([0,0],p.R3);s.m.getIndexOnCircle([0,0],p.rsh)];
    kr = s.m.getIndexOnCircle([0,0],p.R1+p.g/3);
    ks = s.m.getIndexOnCircle([0,0],p.R1+p.g);
    s.m = s.m.makeAG(kr,ks);
    
%     s.m = s.m.getstAMR(find(s.m.tzi(:,end)));
    s.m.showmeshfb
    F(i) = getframe(gcf);
    
    
    s.m = s.m.evalKeFeC('TL3');
    s = s.clearallbcs;
    s = s.setdbc(k0,0);
    for j = 1:length(cc)
        % phase A
        ia = cc(j);
        for k = [1 8 13 20]
            s = s.setExcitation(['c',num2str(k)],ia*Np,'C');
        end
        for k = [2 7 14 19]
            s = s.setExcitation(['c',num2str(k)],-ia*Np,'C');
        end
        s = s.assignEdata;
        s = s.solve(1e-6,10);
        s = s.evalB;
        s = s.evalEtot;
        w(i,j) = s.Etot;
        aa(i,j) = s.evalLF('c1')-s.evalLF('c2')+s.evalLF('c8')-s.evalLF('c7');
%         T(i,j) = s.evalMSTOnAG;
    end
end

close all
hold all
for i=1:length(tt)
    for j = 1:length(cc)
        plot(cc,aa(i,:)*55*188/1000);
%         plot(cc,aaa(i,:)*2*55*188/1000);
%         plot(cc,aa(i,:)*2*55*188/1000,'.','color','k');
    end
end
xlabel('Phase Current [A]');
ylabel('Linkage Flux [wb]')

% close all
% hold all
% teta = linspace(0,45,length(tt));
% for i=1:length(tt)
%     for j = 1:length(cc)
%         plot(teta,T(:,j)*55*1e-9/0.2/(4*pi*1e-7));
% %         torque = spline(teta,w(:,j));
% %         torque.coefs = torque.coefs*diag(3:-1:1,1);
% %         plot(teta,ppval(torque,teta)/10);
%     end
% end