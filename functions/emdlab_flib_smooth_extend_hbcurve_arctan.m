function [h_full,b_full] = emdlab_flib_smooth_extend_hbcurve_arctan(h, b, h_final)
% emdlab_flub_smooth_extend_hbcurve_arctan
% extend b-h curve to h_final using arctangent, smoothly decreasing
% incremental permeability from last slope to mu0.

    mu0 = 4*pi*1e-7;  % vacuum permeability

    % last point
    Hs = h(end);
    Bs = b(end);

    % incremental permeability at last point
    mu_last = (b(end) - b(end-1)) / (h(end) - h(end-1));

    % arctangent scale parameter to control slope decrease
    Ht = linspace(Hs, h_final, 100);  % extension points
    Ht = Ht(2:end);  % skip Hs

    % compute amplitude so that slope decreases from mu_last -> mu0
    % dB/dH = mu0 + Bk*(2/pi)*(1/(1+(H/Hk)^2))*(1/Hk)
    % at H = Hs, slope = mu_last => Bk formula:
    Hk = (h_final - Hs)/2;           % characteristic width, adjustable
    Bk = (mu_last - mu0) * (pi/2) * Hk;  % ensures slope starts at mu_last

    % arctangent extension
    b_ext = Bs + mu0*(Ht - Hs) + Bk*(2/pi)*(atan(Ht/Hk) - atan(Hs/Hk));

    % combine original + extension
    h_full = [h; Ht(:)];
    b_full = [b; b_ext(:)];
end
