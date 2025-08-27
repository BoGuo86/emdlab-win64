function y = emdlab_flib_ppval(p, xq)

y = xq;
x_min = p.breaks(1);
x_max = p.breaks(end);

% internal points
index = (xq >= x_min) & (xq <= x_max);
y(index) = ppval(p, xq(index));

% left points
index = (xq < x_min);
c = p.coefs(1,:);
slope = 3*c(1)*x_min^2 + 2*c(2)*x_min + c(3);
y(index) = ppval(p, x_min) + slope*(xq(index)-x_min);

% right points
index = (xq > x_max);
c = p.coefs(end,:);
slope = 3*c(1)*x_max^2 + 2*c(2)*x_max + c(3);
slope = (ppval(p, p.breaks(end)) - ppval(p, p.breaks(end-1)))/(p.breaks(end)-p.breaks(end-1));
y(index) = ppval(p, x_max) + slope*(xq(index)-x_max);