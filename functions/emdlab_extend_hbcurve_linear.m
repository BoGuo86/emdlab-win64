
function [h,b] = emdlab_extend_hbcurve_linear(h, b, h_final)

mu0 = 4*pi*1e-7*1.1;
x_tmp = linspace(h(end),h_final,1e2)';
y_tmp = b(end) + mu0 * (x_tmp - h(end));

h = [h;x_tmp(2:end)];
b = [b;y_tmp(2:end)];

end