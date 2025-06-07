

clc
clear

t = linspace(0, pi/2, 30);

ip = 20*[cos(t);sin(t)]';
op = 23*[cos(t);sin(t)]';

m = mc_arcAirGap(ip, op)
