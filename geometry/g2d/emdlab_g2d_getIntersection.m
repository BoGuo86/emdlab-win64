function [n,ip,t] = emdlab_g2d_getIntersection(obj1, obj2)

if isa(obj1,'emdlab_g2d_line') && isa(obj2,'emdlab_g2d_line')
    [n,ip,t] = getLineLineIntersection(obj1, obj2);
else
end

end

function [n,ip,t] = getLineLineIntersection(l1, l2)

arguments
    l1 (1,1) emdlab_g2d_line
    l2 (1,1) emdlab_g2d_line
end

d0 = l1.u;
d1 = l2.u;
d1.rotateAroundOrigin(pi/2);
d = l2.p0 - l1.p0;
denom = d0.dot(d1);
d0.rotateAroundOrigin(pi/2);

% default outputs
t = [0,0];
ip = [];

if denom ~= 0

    t(1) = d.dot(d1) / denom;
    t(2) = d.dot(d0) / denom;
    ip = l1.getPoint(t(1));
    n = 1;

else

    tmp = l2.classifyPoint(l1.p0);
    if tmp == emdlab_g2d_location.pl_on
        n = -1;
    else
        n = 0;
    end

end

end