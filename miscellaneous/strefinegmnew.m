
function strefinegmnew()

global mdb gdb

tic
disp('*******************************************************')

pair = [mdb.gm.ConnectivityList(:,[1,2]);
    mdb.gm.ConnectivityList(:,[1,3]);
    mdb.gm.ConnectivityList(:,[2,3])];

[pair,~,ic] = unique(sort(pair,2),'rows');
pmid=(mdb.gm.Points(pair(:,1),:)+mdb.gm.Points(pair(:,2),:))/2;

t1=mdb.gm.ConnectivityList(:,1);
t2=mdb.gm.ConnectivityList(:,2);
t3=mdb.gm.ConnectivityList(:,3);
t12=ic(1:mdb.Ngt)+mdb.Ngp;
t13=ic(mdb.Ngt+1:2*mdb.Ngt)+mdb.Ngp;
t23=ic(2*mdb.Ngt+1:3*mdb.Ngt)+mdb.Ngp;
    
gt=[t1,t12,t13;
    t12,t23,t13;
    t2,t23,t12;
    t3,t13,t23];

gp=[mdb.gm.Points;pmid];
    
% update global triangulation
mdb.gm = triangulation(gt,gp);

% triangle zone index of new mesh
mdb.tzi = repmat(mdb.tzi,1,4);
mdb.tzi = mdb.tzi(:);


toc
tic
% finding index of each boundary points of segment in global mesh
bss = mdb.bs.keys;

for i = 1:numel(bss)
    bs = mdb.bs(bss{i});
    p = gp(bs.pi,:);% points of segment on old mesh
    p = (p(1:end-1,:)+p(2:end,:))/2;
    [~,inp] = ismembertol(p,gp,gdb.geps,'ByRows',true);
    temp = [bs.pi(1:end-1),inp];
    temp = temp';
    bs.pi = [temp(:);bs.pi(end)];
    mdb.bs(bss{i}) = bs;
end

toc
tic
mdb.Ngp = size(mdb.gm.Points,1);
mdb.Ngt = size(mdb.gm.ConnectivityList,1);
mdb.Nge = size(mdb.gm.edges,1);

fprintf('global mesh refined\n')
fprintf('total mesh points = %d\n',mdb.Ngp)
fprintf('total mesh triangles = %d\n',mdb.Ngt)
fprintf('number of mesh zones = %d\n',numel(mdb.mzs))
fprintf('number of boundaries = %d\n',numel(mdb.bs.keys))

toc

end