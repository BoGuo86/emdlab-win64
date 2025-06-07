
tic
clc
clear

global gdb mdb sdb

gdb = gdbc;
mdb = mdbc;
sdb = sdbc;
h = 0.1;

c = [1,0];
newdarcwdkps(c,[2,1],[0,1],0.9*h)
c = [1,1];
newdarcwdkps(c,[2,1],[0,1],0.9*h)
c = [1,-1];
newdarcwdkps(c,[2,1],[0,1],0.9*h)

newcbwd('d',{'A1','A2'},[1 -1])

newcbwd('d1',{'A1','A3'},[1 -1])


[p,t,ss,ipos] = gmesh('d',h);
newmz('rmesh',p,t,ss,ipos,'iron');

[p,t,ss,ipos] = gmesh('d1',h);
newmz('r1mesh',p,t,ss,ipos,'iron');

ggmesh
plotgm