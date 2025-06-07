clc
clear

g = GDBC2D;

g = g.newdlinewdkps([0,0],[2,0],'l1',0.05,'l2',0.1);
g = g.newdlinewdkps([2,0],[2,2],'l1',0.1,'l2',0.2);
g = g.newdlinewdkps([0,2],[2,2],'l1',0.1,'l2',0.2);
g = g.newdlinewdkps([0,0],[0,2],'l1',0.05,'l2',0.1);
g = g.newcb('s','L1',1,'L2',1,'L3',-1,'L4',-1);

g = g.newdDM('s','s');
% g.plotwf

m = MDBCT