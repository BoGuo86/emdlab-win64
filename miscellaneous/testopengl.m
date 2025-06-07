

% p = 3*rand(100,2)-1;
% t = delaunay(p);

p = m.mzs.stator.nodes./max(max(m.mzs.stator.nodes));
t = m.mzs.stator.cl;

oglTest(t,p);

clear mex
