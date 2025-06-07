a = m.mts.steel_1008.vB;
x = linspace(0,5,50);
x_tmp = x;
x_tmp(x>a.breaks(end)) = a.breaks(end);
plot(x,ppval(a, x_tmp))