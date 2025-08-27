function [df_di, df_dj, df_dk] = emdlab_flib_gradient3(f, di, dj, dk)
% Computes the gradient of a 3D scalar field f(i,j,k)
% Inputs:
%   f  - 3D array (n1 × n2 × n3)
%   di - spacing in i (x-direction)
%   dj - spacing in j (y-direction)
%   dk - spacing in k (z-direction)
% Outputs:
%   df_di - partial derivative ∂f/∂i
%   df_dj - partial derivative ∂f/∂j
%   df_dk - partial derivative ∂f/∂k

[n1, n2, n3] = size(f);

df_di = zeros(size(f));  % ∂f/∂i
df_dj = zeros(size(f));  % ∂f/∂j
df_dk = zeros(size(f));  % ∂f/∂k

for i = 1:n1
    for j = 1:n2
        for k = 1:n3

            % ∂f/∂i (x-direction)
            if i == 1
                df_di(i,j,k) = (f(i+1,j,k) - f(i,j,k)) / di;
            elseif i == n1
                df_di(i,j,k) = (f(i,j,k) - f(i-1,j,k)) / di;
            else
                df_di(i,j,k) = (f(i+1,j,k) - f(i-1,j,k)) / (2*di);
            end

            % ∂f/∂j (y-direction)
            if j == 1
                df_dj(i,j,k) = (f(i,j+1,k) - f(i,j,k)) / dj;
            elseif j == n2
                df_dj(i,j,k) = (f(i,j,k) - f(i,j-1,k)) / dj;
            else
                df_dj(i,j,k) = (f(i,j+1,k) - f(i,j-1,k)) / (2*dj);
            end

            % ∂f/∂k (z-direction)
            if k == 1
                df_dk(i,j,k) = (f(i,j,k+1) - f(i,j,k)) / dk;
            elseif k == n3
                df_dk(i,j,k) = (f(i,j,k) - f(i,j,k-1)) / dk;
            else
                df_dk(i,j,k) = (f(i,j,k+1) - f(i,j,k-1)) / (2*dk);
            end

        end
    end
end

end
