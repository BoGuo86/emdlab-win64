clc
clear

g = GDBC2D;

g = g.newdlinewdkps([0,0],[1,0],'nnodes',8,'alpha',0.8);
g = g.newdarcpl('L1','p2','f',1,90,-1);
g = g.newdlinewdkps([2,1],[2,2],'nnodes',8,'alpha',1.2);
g = g.newdlinewdkps([2,2],[0,2],'nnodes',8);
g = g.newdlinewdkps([0,2],[0,0],'nnodes',8);

g = g.newcb('b','L1',1,'A1',1,'L2',1,'L3',1,'L4',1);

g = g.newdDM('s','b');
% g = g.newdd('s',0.2,'b');

% g = g.newpolygonald('s',inf,[0,0;1,0;1,1;0,1]);
% 
% g.plotinitmesh
% 
m = MDBCT(g);


% m = m.ggmesh;
% m = m.strefine;
% m.showmesh

% g = g.newdlinewdkps([0,0],[1,1]);


