function [B_ext,nu_ext] = emdlab_extend_linear_to_air(B,nu)

mu0 = 4*pi*1e-7;
nu0 = 1/mu0;

% derivative at last point
df = fnder(pchip(B,nu));
slope = ppval(df, B(end));

% distance needed to reach air reluctivity
DeltaB = (nu0 - nu(end)) / slope;

% extend moderately
B_ext = [B;
         B(end) + DeltaB;
         B(end) + 10*DeltaB];

nu_ext = [nu;
          nu0;
          nu0];

end
