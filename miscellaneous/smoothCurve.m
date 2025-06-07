function y = smoothCurve(x)
y = x;
for i = 1:size(x,1)
    y(i,:) = smooth(x(i,:));
end
end

function y = smooth(x)
n = 1:length(x);
p = pchip(n,x);
xx = n+0.5;
xx = xx(1:end-1);
yy = ppval(p, xx);
p = pchip(xx,yy);
y = ppval(p, n);
end