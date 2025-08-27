function [x, y] = emdlab_g2d_getLineLineIntersection(p1, u1, p2, u2)

A = [u1(1),-u2(1);u1(2),-u2(2)];
b = [p2(1)-p1(1);p2(2)-p1(2)];
u = A\b;
x = p1(1) + u(1) * u1(1);
y = p1(2) + u(1) * u1(2);

end