
[x,y] = meshgrid(linspace(-1.5,-0.5,3));
m = TMDBC;
m = m.addmz('s1',TMZPC(delaunay([x(:),y(:)]),[x(:),y(:)]));
m = m.cshiftmz('s2','s1',[2,0]);
m = m.cshiftmz('s3','s1',[2,2]);
m = m.cshiftmz('s4','s1',[0,2]);
m = m.ggmesh;
m = m.makeRegion([-3,-3],[3,3],8,8);
m = m.strefine;
m = m.strefine;
m = m.strefine;
m.showmzs
% v = 2*[5,0;10,0;10,9.5;9.5,10;9,10;9,1;5,1];
% f = [1:size(v,1);[2:size(v,1),1]]';
% 
% t = delaunayTriangulation(v,f);
% t = t.ConnectivityList(ext_dpoly2d(t.incenter,f,v)<0,:);
% m = TMZPC(t,v);
% m = m.getRotateY(2*pi,60);
% m = m.setdata;
% m.showwf
% m = TTMDBC();
% m = m.addmz('m1',getRotateY(TMZPC(y2,y1),pi/2,20));
% m = m.addmz('m2',getRotateX(TMZPC(y2,y1),-pi/2,20));
% m = m.ggmesh;
% m.showwf


