
function [B,v] = extend_exponential(B, v, B_final)

slope = (v(end)-v(end-1))/(B(end)-B(end-1));
tau = (B_final-v(end))/slope;

x_tmp = linspace(B(end),1e2*B(end),1e2)';
y_tmp = B_final + (v(end)-B_final) * exp(-(x_tmp-B(end))/tau);

B = [B;x_tmp(2:end)];
v = [v;y_tmp(2:end)];

end