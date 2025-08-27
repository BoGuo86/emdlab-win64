function df_dx = emdlab_flib_gradient1(f, dx)
% Calculate the gradient (derivative) of a 1D array f with respect to x
% Inputs:
%   f  - 1D array of function values
%   dx - spacing between points
% Output:
%   df_dx - numerical derivative of f with respect to x

n = length(f);
df_dx = zeros(size(f));

for i = 1:n
    if i == 1
        df_dx(i) = (f(i+1) - f(i)) / dx;  % forward difference
    elseif i == n
        df_dx(i) = (f(i) - f(i-1)) / dx;  % backward difference
    else
        df_dx(i) = (f(i+1) - f(i-1)) / (2*dx);  % central difference
    end
end

end
