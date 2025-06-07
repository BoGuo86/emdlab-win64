
function ggmeshnew()

global mdb gdb

% generation of global mesh
tic
disp('*******************************************************')
% total points in global mesh
gp = zeros(mdb.Np,2);

% total triangles in global mesh
gt = zeros(mdb.Nt,3);

% zone index of each triangle
mdb.tzi = zeros(mdb.Nt,1);

% getting all mesh zones data
mzs = mdb.mzs.values;
temp = 0;
for i = 1:numel(mzs)
    % number of ith zone points
    Np = size(mzs{i}.p,1);
    gp(temp+1:temp+Np,:)  = mzs{i}.p;
    temp = temp + Np;
end

% unique global points
[gp,~,ic] = uniquetol(gp,gdb.geps,'ByRows',true);

% forming total triangles of global mesh
temp = 0;
tempi = 0;
for i = 1:numel(mzs)
    % number of ith zone points
    Np = size(mzs{i}.p,1);
    mzs{i}.gpi = ic(tempi+1:Np+tempi);
    
    % getting indices of ith zone boundary points in global mesh
    for j = 1:numel(mzs{i}.bs)      
        bs = M_bspc;
        bs.pi = mzs{i}.gpi(mzs{i}.ipos{j});
        mdb.bs(mzs{i}.bs{j}) = bs;
    end
    
    % number of ith zone triangles
    Nt = size(mzs{i}.t,1);
    
    % adding new triangles to global triangle list
    gt(temp+1:temp+Nt,:)  = mzs{i}.gpi(mzs{i}.t);
    
    % setting triangles zone index
    mdb.tzi(temp+1:temp+Nt,1) = mzs{i}.zi;
    temp = temp + Nt;
    tempi = tempi + Np;
end

% needed data
mdb.mzs = containers.Map(mdb.mzs.keys,mzs);
mdb.gm = triangulation(gt,gp);
mdb.Ngp = size(mdb.gm.Points,1);
mdb.Ngt = size(mdb.gm.ConnectivityList,1);
mdb.Nge = size(mdb.gm.edges,1);

% specefing mesh zones
tzi = zeros(mdb.Ngt,mdb.Nmzs);
for i = 1:mdb.Nmzs
    tzi(:,i) = mdb.tzi(:,1) == i; 
end
mdb.tzi = logical(tzi);
    

fprintf('global mesh created\n')
fprintf('total mesh points = %d\n',mdb.Ngp)
fprintf('total mesh triangles = %d\n',mdb.Ngt)
fprintf('number of mesh zones = %d\n',numel(mdb.mzs))
fprintf('number of boundaries = %d\n',numel(mdb.bs.keys))
toc

end


