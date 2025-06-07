function [e,p] = tri2pelements(t,p,z)

if iscolumn(z)
    z = z';
end

Np = size(p,1);
Nz = length(z);
p = repmat(p,length(z),1);
z = repmat(z,Np,1);
p = [p,z(:)];

ip1 = t(:,1);
ip2 = t(:,2);
ip3 = t(:,3);
ip4 = ip1+Np;
ip5 = ip2+Np;
ip6 = ip3+Np;

e1 = [ip1,ip2,ip3,ip4,ip5,ip6];

e = zeros([],6);

for i = 1:Nz-1
    e = [e;e1+(i-1)*Np];
end

end