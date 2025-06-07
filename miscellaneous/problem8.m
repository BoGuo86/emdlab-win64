
x = [-3,-2,-1,0,1,2,3];
y = [0,0,0,1,0,0,0];

xq = -3:.01:3;
p = pchip(x,y,xq);
s = spline(x,y,xq);
plot(x,y,'o',xq,p,'-',xq,s,'-.')
legend('Sample Points','pchip','spline','Location','SouthEast')
