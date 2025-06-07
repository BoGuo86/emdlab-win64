
function [x,y] = extend_exponential(x, y, y_final)

% p = spline(x,y);
% c = p.coefs(end,:);
% x_tmp = x(end) - x(end-1);
% slope = 3*c(4)*x_tmp^2 + 2*c(3)*x_tmp + c(2);

slope = (y(end)-y(end-1))/(x(end)-x(end-1));
tau = (y_final-y(end))/slope;

x_tmp = linspace(x(end),1e2*x(end),1e2)';
y_tmp = y_final + (y(end)-y_final) * exp(-(x_tmp-x(end))/tau);

x = [x;x_tmp(2:end)];
y = [y;y_tmp(2:end)];

end