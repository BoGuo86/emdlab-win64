p = [-1,-0.5;1,-1.5;1,1.5;-1,0.5];
t = [1,2,3;1,3,4];

m = TMZPC(t,p);

m = m.strefine;
m = m.moveNodes;
m = m.strefine;
m = m.moveNodes;
m = m.strefine;
m = m.moveNodes;
m = m.strefine;
m = m.moveNodes;

m.showm