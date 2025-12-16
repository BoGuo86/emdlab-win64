function [B_ext, v_ext] = emdlab_flib_smooth_extend_bvcurve_arctan(B, v, v_final)
% emdlab_flib_smooth_extend_bvcurve_arctan
% extend B and reluctivity v to v_final using smooth arctangent decay
%
% inputs:
%   B       : original B vector (t)
%   v       : original reluctivity vector (1/mu or v=B/H)
%   v_final : asymptotic v value for high B
%
% outputs:
%   B_ext   : extended B vector
%   v_ext   : extended v vector

    % last point
    Bs = B(end);
    vs = v(end);

    % last slope
    slope_last = (v(end) - v(end-1)) / (B(end) - B(end-1));

    % characteristic B width for arctan transition
    Bk = (B(end)*0.5);  % you can tune this fraction

    % compute amplitude to match slope at B=end
    A = (slope_last - 0) * (pi/2) * Bk;  % slope decays from slope_last to near 0

    % generate extension points
    B_tail = linspace(B(end), 1e2*B(end), 100)'; 
    B_tail = B_tail(2:end); % skip first point

    % arctangent decay
    v_tail = v_final + A * (atan((B_tail - Bs)/Bk) - atan(0));

    % combine original + extension
    B_ext = [B; B_tail];
    v_ext = [v; v_tail];

end
