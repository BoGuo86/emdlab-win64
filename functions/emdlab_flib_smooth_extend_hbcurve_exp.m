function [h_ext, b_ext] = emdlab_flib_smooth_extend_hbcurve_exp(h, b, h_final)
% emdlab_flib_smooth_extend_hbcurve_exp
% smooth maxwell-compatible extension of b-h curve
%
% inputs:
%   h       : original h vector (a/m), monotonic increasing
%   b       : original b vector (t)
%   h_final : final h value for extension (a/m)
%
% outputs:
%   h_ext   : extended h vector
%   b_ext   : extended b vector

    % constants
    mu0 = 4*pi*1e-7;

    % ensure column vectors
    h = h(:);
    b = b(:);

    % last data point (saturation entry)
    Hsat = h(end);
    Bsat = b(end);

    % incremental permeability at last segment
    mu_last = (b(end) - b(end-1)) / (h(end) - h(end-1));

    % smoothing length (maxwell-like)
    DeltaH = 100 * Hsat;

    % interpolation inside measured range
    pp = pchip(h, b);

    % generate extended h grid
    if h_final <= Hsat
        h_ext = h;
        b_ext = b;
        return
    end

    dh = mean(diff(h));
    h_tail = (Hsat+dh : dh : h_final).';
    h_ext  = [h; h_tail];

    % allocate b
    b_ext = zeros(size(h_ext));

    % inside data range
    idx_in = h_ext <= Hsat;
    b_ext(idx_in) = ppval(pp, h_ext(idx_in));

    % smooth asymptotic extension
    idx_out = h_ext > Hsat;
    Ht = h_ext(idx_out) - Hsat;
    b_ext(idx_out) = Bsat + mu0 * Ht + (mu_last - mu0) * DeltaH .* (1 - exp(-Ht / DeltaH));

end
