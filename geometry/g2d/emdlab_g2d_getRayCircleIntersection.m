function [x,y] = emdlab_g2d_getRayCircleIntersection(p, u, c, r)

u = u/norm(u);
err = @(t) abs(norm(p+t*u-c)-r);
options = optimset('TolX',1e-6);
t = fminbnd(err,0,r,options);
p0 = p+t*u;
x = p0(1);
y = p0(2);

end