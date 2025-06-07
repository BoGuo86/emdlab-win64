function [np,nq] = strefineQ(p,q)

% standard refinement
% first edge of all triangles
e1 = q(:,[1,2]);
% second edge of all triangles
e2 = q(:,[2,3]);
% third edge of all triangles
e3 = q(:,[3,4]);
% forth edge of all triangles
e4 = q(:,[4,1]);
% total edges
e = [e1;e2;e3;e4];

% middle point of each edge
pmid = (p(e(:,1),:)+p(e(:,2),:))/2;

% center points of each quadrilateral
pc = (p(q(:,1),:)+p(q(:,2),:)+p(q(:,3),:)+p(q(:,4),:))/4;

% making unique new points
[pmid,~,ic] = uniquetol(pmid,1e-3,'ByRows',true);

% total unique new points
np = [p;pmid;pc];

% number of old points
Np = size(p,1);

% number of old triangles
Nq = size(q,1);

% index of point 1 of each triangle
ip1 = q(:,1);

% index of point 2 of each triangle
ip2 = q(:,2);

% index of point 3 of each triangle
ip3 = q(:,3);

% index of point 4 of each triangle
ip4 = q(:,4);

% index of point 5 of each triangle
ip5 = ic(1:Nq)+Np;

% index of point 6 of each triangle
ip6 = ic(Nq+1:2*Nq)+Np;

% index of point 7 of each triangle
ip7 = ic(2*Nq+1:3*Nq)+Np;

% index of point 8 of each triangle
ip8 = ic(3*Nq+1:4*Nq)+Np;

% index of point 9 of each triangle
ip9 = (1:Nq)'+Np+size(pmid,1);

% new triangles
nq = [ip1,ip5,ip9,ip8
    ip5,ip2,ip6,ip9
    ip9,ip6,ip3,ip7
    ip8,ip9,ip7,ip4];

end
