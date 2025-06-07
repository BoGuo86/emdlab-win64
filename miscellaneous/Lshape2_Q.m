
p = [-1,0
    1,0
    2,0
    2,1
    2,3
    1,3
    1,1
    0,1];

q = [1 2 7 8
    2 3 4 7
    4 5 6 7];

plotQmesh(p,q)

%%
[p,q] = strefineQ(p,q);
[p,q] = strefineQ(p,q);
[p,q] = strefineQ(p,q);
plotQmesh(p,q)