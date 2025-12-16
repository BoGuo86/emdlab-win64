function [B_ext, v_ext] = emdlab_flib_smooth_extend_bvcurve_exp(B, v, B_final)
% extend_exponential
% extend B and reluctivity v to B_final using smooth exponential decay
%
% inputs:
%   B       : original B vector (T)
%   v       : original reluctivity vector (1/mu or v=B/H)
%   B_final : asymptotic B value for v
%
% outputs:
%   B_ext   : extended B vector
%   v_ext   : extended v vector

    % compute slope from last segment
    slope = (v(end) - v(end-1)) / (B(end) - B(end-1));

    % characteristic decay constant
    tau = (B_final - v(end)) / slope;

    % create extension points
    B_tail = linspace(B(end), 1e2*B(end), 100)';  % 100 points beyond last B
    B_tail = B_tail(2:end); % skip first point (already included)

    % exponential decay for v
    v_tail = B_final + (v(end) - B_final) * exp(-(B_tail - B(end))/tau);

    % combine original + extension
    B_ext = [B; B_tail];
    v_ext = [v; v_tail];

end
