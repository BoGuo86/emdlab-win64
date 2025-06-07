clc
clear

% main dim
rsh = 8.000000;
R0 = 30.400000;
R1 = 40.900000;
R2 = 65.800000;
R3 = 74.900000;
betar = 19 * pi/180;
betas = 16 * pi/180;
g = 0.3; 

gammar = asin(R1*sin(betar/2)/R0);
gammas = asin((R1+g)*sin(betas/2)/R2);

thetag = 1;

Kh = [rsh^2,rsh,1;
    (R1+g/2)^2,R1+g/2,1;
    R3^2,R3,1];
ch = Kh\[2;1;2];

fh = @(p) ch(1)*(sqrt(sum(p.^2,2))).^2+ch(2)*sqrt(sum(p.^2,2))+ch(3);

% key points

rkp = zeros(12,2);
rkp(1,:) = rsh*[cos(-pi/4),sin(-pi/4)];
rkp(2,:) = R0*[cos(-pi/4),sin(-pi/4)];
rkp(11,:) = R0*[cos(pi/4),sin(pi/4)];
rkp(12,:) = rsh*[cos(pi/4),sin(pi/4)];

rkp(3,:) = R0*[cos(gammar),sin(-gammar)];
rkp(4,:) = R1*[cos(betar/2),sin(-betar/2)];
rkp(5,:) = R1*[cos(betar/2),sin(betar/2)];
rkp(6,:) = R0*[cos(gammar),sin(gammar)];

rkp(3:6,:) = protate(rkp(3:6,:),-pi/8);
rkp(7:10,:) = protate(rkp(3:6,:),pi/4);

% plot(rkp(:,1),rkp(:,2),'*');axis off square

% mesh max length
hr = 1;
hs = 1;
hg = 1;

s.rs1 = getlinepoints(rkp(1,:),rkp(2,:),0.9*hr);
s.rs2 = getarcpoints([0,0],rkp(2,:),rkp(3,:),0.9*min(hr,hg),5);
s.rs3 = getlinepoints(rkp(3,:),rkp(4,:),0.9*min(hr,hg));
s.rs4 = getarcpoints([0,0],rkp(4,:),rkp(5,:),0.9*min(hr,hg),thetag);
s.rs5 = getlinepoints(rkp(5,:),rkp(6,:),0.9*min(hr,hg));
s.rs6 = getarcpoints([0,0],rkp(6,:),rkp(7,:),0.9*min(hr,hg),5);
s.rs7 = getlinepoints(rkp(7,:),rkp(8,:),0.9*min(hr,hg));
s.rs8 = getarcpoints([0,0],rkp(8,:),rkp(9,:),0.9*min(hr,hg),thetag);
s.rs9 = getlinepoints(rkp(9,:),rkp(10,:),0.9*min(hr,hg));
s.rs10 = getarcpoints([0,0],rkp(10,:),rkp(11,:),0.9*min(hr,hg),5);
s.rs11 = getlinepoints(rkp(11,:),rkp(12,:),0.9*hr);
s.rs12 = getarcpoints([0,0],rkp(12,:),rkp(1,:),0.9*2*hr,15);


rpoly = [rkp(1,:);s.rs1;rkp(2,:);s.rs2;rkp(3,:);s.rs3;rkp(4,:);s.rs4;...
    rkp(5,:);s.rs5;rkp(6,:);s.rs6;rkp(7,:);s.rs7;rkp(8,:);s.rs8;...
    rkp(9,:);s.rs9;rkp(10,:);s.rs10;rkp(11,:);s.rs11;rkp(12,:);s.rs12;...
    rkp(1,:)];

gkp1 = (R1+g/2)*[cos(pi/4),sin(pi/4)];
gkp2 = (R1+g/2)*[cos(-pi/4),sin(-pi/4)];
s.gs = getarcpoints([0,0],gkp1,gkp2,hg,thetag);
gsl = getlinepoints(gkp2,rkp(2,:),0.9*hg);
gsu = getlinepoints(rkp(11,:),gkp1,0.9*hg);

g1poly = [gkp1;s.gs;gkp2;gsl;rkp(2,:);...
    s.rs2;rkp(3,:);s.rs3;rkp(4,:);s.rs4;...
    rkp(5,:);s.rs5;rkp(6,:);s.rs6;rkp(7,:);s.rs7;rkp(8,:);s.rs8;...
    rkp(9,:);s.rs9;rkp(10,:);s.rs10;rkp(11,:);gsu;gkp1];

% plot(g1poly(:,1),g1poly(:,2))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
skp = zeros(22,2);

skp(1,:) = R3*[cos(-pi/4),sin(-pi/4)];
skp(2,:) = R2*[cos(-pi/4),sin(-pi/4)];
skp(18,:) = R3*[cos(pi/4),sin(pi/4)];
skp(17,:) = R2*[cos(pi/4),sin(pi/4)];
skp(19,:) = (R1+g)*[cos(-pi/4),sin(-pi/4)];
skp(22,:) = (R1+g)*[cos(pi/4),sin(pi/4)];
skp(7,:) = protate(skp(2,:),pi/6);
skp(12,:) = protate(skp(7,:),pi/6);
skp(20,:) = protate(skp(19,:),pi/6);
skp(21,:) = protate(skp(20,:),pi/6);

skp(10,:) = (R1+g)*[cos(betas/2),sin(betas/2)];
skp(9,:) = [skp(10,1),-skp(10,2)];
skp(11,:) = R2*[cos(gammas),sin(gammas)];
skp(8,:) = [skp(11,1),-skp(11,2)];

skp(3:6,:) = protate(skp(8:11,:),-pi/6);
skp(13:16,:) = protate(skp(8:11,:),pi/6);

% plot(skp(:,1),skp(:,2),'*');axis off square

ss1 = getlinepoints(skp(1,:),skp(2,:),0.9*1.5*hs);
ss2 = getarcpoints([0,0],skp(2,:),skp(3,:),0.9*min(hs,hg),5);
ss3 = getlinepoints(skp(3,:),skp(4,:),0.9*min(hs,hg));
ss4 = getarcpoints([0,0],skp(4,:),skp(5,:),0.9*min(hs,hg),thetag);
ss5 = getlinepoints(skp(5,:),skp(6,:),0.9*min(hs,hg));
ss6 = getarcpoints([0,0],skp(6,:),skp(7,:),0.9*min(hs,hg),5);
ss7 = getarcpoints([0,0],skp(7,:),skp(8,:),0.9*min(hs,hg),5);
ss8 = getlinepoints(skp(8,:),skp(9,:),0.9*min(hs,hg));
ss9 = getarcpoints([0,0],skp(9,:),skp(10,:),0.9*min(hs,hg),thetag);
ss10 = getlinepoints(skp(10,:),skp(11,:),0.9*min(hs,hg));
ss11 = getarcpoints([0,0],skp(11,:),skp(12,:),0.9*min(hs,hg),5);
ss12 = getarcpoints([0,0],skp(12,:),skp(13,:),0.9*min(hs,hg),5);
ss13 = getlinepoints(skp(13,:),skp(14,:),0.9*min(hs,hg));
ss14 = getarcpoints([0,0],skp(14,:),skp(15,:),0.9*min(hs,hg),thetag);
ss15 = getlinepoints(skp(15,:),skp(16,:),0.9*min(hs,hg));
ss16 = getarcpoints([0,0],skp(16,:),skp(17,:),0.9*min(hs,hg),5);
ss17 = getlinepoints(skp(17,:),skp(18,:),0.9*1.5*hs);
ss18 = getarcpoints([0,0],skp(18,:),skp(1,:),0.9*2*hs,10);

spoly = [skp(1,:);ss1;
    skp(2,:);ss2;
    skp(3,:);ss3;
    skp(4,:);ss4;
    skp(5,:);ss5;
    skp(6,:);ss6;
    skp(7,:);ss7;
    skp(8,:);ss8;
    skp(9,:);ss9;
    skp(10,:);ss10;
    skp(11,:);ss11;
    skp(12,:);ss12;
    skp(13,:);ss13;
    skp(14,:);ss14;
    skp(15,:);ss15;
    skp(16,:);ss16;
    skp(17,:);ss17;
    skp(18,:);ss18;skp(1,:)];

% plot(spoly(:,1),spoly(:,2),'*');axis off square


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ss19 = getlinepoints(skp(19,:),skp(2,:),0.9*min(hs,hg));
ss20 = getlinepoints(skp(20,:),skp(7,:),0.9*min(hs,hg));
ss21 = getlinepoints(skp(21,:),skp(12,:),0.9*min(hs,hg));
ss22 = getlinepoints(skp(22,:),skp(17,:),0.9*min(hs,hg));

ss23 = getarcpoints([0,0],skp(19,:),skp(4,:),0.9*min(hs,hg),thetag);
ss24 = getarcpoints([0,0],skp(5,:),skp(20,:),0.9*min(hs,hg),thetag);
ss25 = getarcpoints([0,0],skp(20,:),skp(9,:),0.9*min(hs,hg),thetag);
ss26 = getarcpoints([0,0],skp(10,:),skp(21,:),0.9*min(hs,hg),thetag);
ss27 = getarcpoints([0,0],skp(21,:),skp(14,:),0.9*min(hs,hg),thetag);
ss28 = getarcpoints([0,0],skp(15,:),skp(22,:),0.9*min(hs,hg),thetag);

%%%%%%%%%
g2poly = [gkp1;s.gs;gkp2;
    skp(19,:);ss23;
    skp(4,:);ss4;
    skp(5,:);ss24;
    skp(20,:);ss25;
    skp(9,:);ss9;
    skp(10,:);ss26;
    skp(21,:);ss27;
    skp(14,:);ss14;
    skp(15,:);ss28;
    skp(22,:);gkp1];

% plot(g2poly(:,1),g2poly(:,2))



%%%%%%%%%%%%%

c1a1poly = [skp(10,:);ss10;
    skp(11,:);ss11;
    skp(12,:);preverse(ss21);
    skp(21,:);preverse(ss26);skp(10,:)];
% plot(c1a1poly(:,1),c1a1poly(:,2))

c1a2poly = [skp(7,:);ss7;
    skp(8,:);ss8;
    skp(9,:);preverse(ss25);
    skp(20,:);ss20;skp(7,:)];

c2a1poly = [skp(21,:);ss21;
    skp(12,:);ss12;
    skp(13,:);ss13;
    skp(14,:);preverse(ss27);skp(21,:)];

c2a2poly = [skp(15,:);ss15;
    skp(16,:);ss16;
    skp(17,:);preverse(ss22);
    skp(22,:);preverse(ss28);skp(15,:)];

c3a1poly = [skp(19,:);ss19;
    skp(2,:);ss2;
    skp(3,:);ss3;
    skp(4,:);preverse(ss23);skp(19,:)];

c3a2poly = [skp(5,:);ss5;
    skp(6,:);ss6;
    skp(7,:);preverse(ss20);
    skp(20,:);preverse(ss24);skp(5,:)];

per1 = [rkp(1,:);s.rs1;rkp(2,:);preverse(gsl);gkp2;skp(19,:);preverse(ss19);...
    skp(2,:);preverse(ss1);skp(1,:)];

per2 = [rkp(12,:);preverse(s.rs11);rkp(11,:);gsu;gkp1;skp(22,:);ss22;...
    skp(17,:);ss17;skp(18,:)];

bound = [per1;preverse(ss18);preverse(per2);s.rs12];

pfix = [s.gs;
    rkp;skp;gkp1;gkp2;
    ss1;
    ss2;
    ss3;
    ss4;
    ss5;
    ss6;
    ss7;
    ss7;
    ss9;
    ss10;
    ss12;
    ss13;
    ss14;
    ss15;
    ss16;
    ss17;
    ss18;
    ss19;
    ss20;
    ss21;
    ss22;
    ss23;
    ss24;
    ss25;
    ss26;
    ss27;
    ss28;
    gsu;
    gsl;
    s.rs1;
    s.rs2;
    s.rs3;
    s.rs4;
    s.rs5;
    s.rs6;
    s.rs7;
    s.rs8;
    s.rs9;
    s.rs10;
    s.rs11;
    s.rs12;];



[p,t] = distmesh2d1(@dpoly,fh,hg,[min(pfix);max(pfix)],pfix,bound);

close all
% e1 = 1:size(s.gs,1);
% e2 = [2:size(s.gs,1),1];
% dt = delaunayTriangulation(p,[e1',e2']);
% triplot(dt,'r');axis off equal

