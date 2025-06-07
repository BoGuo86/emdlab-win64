
function y = divideVector(x)

y = (x(1:end-1) + x(2:end))/2;

y = [x(1:end-1);y];
y = y(:);
y = [y', x(end)];

end

