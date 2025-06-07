function [h, b] = smoothHBCurve(HBData, degree)

if nargin < 2
    degree = 2;
elseif nargin > 2
    error('Too many input arguments');
end

if size(HBData, 2) ~= 2
    error('HB data must be a [Nx2] matrix.');
end

% getting h and b
h = HBData(:, 1);
b = HBData(:, 2);

% normalizing h and b
hmax = h(end);
bmax = b(end);
h = [0;h] / hmax;
b = [0;b] / bmax;

% interpolation according to desired degree
% A: regression matrix
% c: coefficient vector
switch degree
       
    case 1
        
        A = [h, -b .* h];
        c = (A' * A) \ (A' * b);
        b = (c(1) * h) ./ (1 + c(3) * h);
        
    case 12
        
        A = [h, h.^2, -b .* h, -b .* (h.^2)];
        c = (A' * A) \ (A' * b);
        b = (c(1) * h + c(2) * h.^2) ./ (1 + c(3) * h + c(4) * h.^2);
        
    case 135
        
        A = [h, h.^3, h.^5, -b .* h, -b .* (h.^3), -b .* (h.^5)];
        c = (A' * A) \ (A' * b);
        b = (c(1) * h + c(2) * h.^3 + c(3) * h.^5) ./ (1 + c(4) * h + c(5) * h.^3 + c(6) * h.^5);
        
    otherwise
end

% rescaling
b = b(2:end)*bmax;
h = h(2:end)*hmax;

end
