function [h_full, b_full] = emdlab_flib_smooth_extend_hbcurve_polyc(h, b, h_final, n)
% emdlab_flib_smooth_extend_hbcurve_polyc
% smooth polynomial-based extension of B-H curve
% decreasing incremental permeability from last slope to mu0
%
% inputs:
%   h       : original h vector (a/m)
%   b       : original b vector (t)
%   h_final : final h value for extension (a/m)
%   n       : polynomial exponent for decay (default 2)
%
% outputs:
%   h_full  : extended h vector
%   b_full  : extended b vector

    if nargin < 4
        n = 2; % default exponent
    end

    mu0 = 4*pi*1e-7; % vacuum permeability

    % last point
    Hs = h(end);
    Bs = b(end);

    % incremental permeability at last segment
    mu_last = (b(end) - b(end-1)) / (h(end) - h(end-1));

    % generate extended H vector
    Ht = linspace(Hs, h_final, 100);
    Ht = Ht(2:end); % skip Hs

    % polynomial decay function (normalized)
    f = (1 - (Ht - Hs)/(h_final - Hs)).^n;

    % incremental permeability along extension
    mu_ext = mu0 + (mu_last - mu0) .* f;

    % integrate mu(H) to get B(H)
    b_ext = Bs + cumtrapz(Ht, mu_ext);

    % combine original + extension
    h_full = [h; Ht(:)];
    b_full = [b; b_ext(:)];
end
