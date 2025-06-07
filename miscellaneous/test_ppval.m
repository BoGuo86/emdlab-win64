
x = [0,2,3,4];
y = [0,2.5,3,4];

plot(x,y);hold on

sp = spline(x,y);
pc = pchip(x,y);

x = linspace(x(1),x(end),100);
plot(x,ppval(sp,x));
plot(x,ppval(pc,x));