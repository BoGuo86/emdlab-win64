function [B,nu] = emdlab_extend_exponential1(B, nu)

mu0 = 4*pi*1e-7;
nu0 = 1/mu0;

f = pchip(B, nu);
f.coefs = f.coefs * diag(3:-1:1, 1);

slope = ppval(f,B(end));
Delta_B = (nu0-nu(end))/slope;

B = [B;B(end)+Delta_B;B(end)+1e3*Delta_B];
nu = [nu;nu0;nu0];

end