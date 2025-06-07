clc
clear all

global gdb mdb

gdb = G_gdbc;
mdb = M_mdbc;

h = 0.25;

newdlinewdkps([0,0],[1,0],0.995*h)
newdarccppwdkps([0,0],[1,0],[0,1],1,0.995*h)
newdlinewdkps([0,1],[0,0],0.995*h)
newcbwd('p',{'L1','A1','L2'},[1 1 1])

[p,t,ss,ipos] = gmesh('p',h);
newmz('pm',p,t,ss,ipos);

ggmeshnew

plotgm

e = mdb.gm.edges;
e = sort(e,2);

fbe = mdb.gm.freeBoundary;
fbe = sort(fbe,2);

e = setdiff(e,fbe,'rows');

t = mdb.gm.ConnectivityList;
p = mdb.gm.Points;
Nt = size(mdb.gm.ConnectivityList,1);

if rem(Nt,2) == 0
    q = zeros(Nt/2,4);
    for counter = 1:Nt/2-1
        ea = mdb.gm.edgeAttachments(e(1,:));
        tt = t(ea{1},:);
        pp = tt(:);
        pp = (unique(pp))';
        pt = p(pp,:);
        qi = convhull(pt);
        q(counter,:) = pp(qi(1:end-1));
        ee = [tt(:,[1,2]);tt(:,[2,3]);tt(:,[3,1])];
        ee = sort(ee,2);
        ee = unique(ee,'rows');
        e = setdiff(e,ee,'rows');
    end
end
q = q(1:end-1,:);
[p,q] = strefineQ(p,q);
[p,q] = strefineQ(p,q);
plotQmesh(p,q)