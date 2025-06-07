function [e,p] = tri2ttrotate(t,p,Angle,Nz)


Np = size(p,1);
p1 = [p,zeros(Np,1)];
p = zeros([],3);

for i = 1:Nz
    p = [p;protatey(p1,(i-1)*Angle/Nz)];
end

ip1 = t(:,1);
ip2 = t(:,2);
ip3 = t(:,3);
ip4 = ip1+Np;
ip5 = ip2+Np;
ip6 = ip3+Np;

e1 = [ip4,ip1,ip3,ip2
    ip6,ip4,ip3,ip2
    ip6,ip5,ip4,ip2];
e = zeros([],4);

for i = 1:Nz-1
    e = [e;e1+(i-1)*Np];
end


end