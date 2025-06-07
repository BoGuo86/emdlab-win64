
clc
clear

t = linspace(0,pi/2,30)';

ips = [cos(t), sin(t)];
ops = 1.2*[cos(t), sin(t)];

ips = circshift(ips,5);

% [a,b]=gts_sortPointsCCW(ips,[0,0])
m = mc_arcAirGap(ips, ops);
m.rotateInner(pi/12)

m.m.showm
