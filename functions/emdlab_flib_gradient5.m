function [df_di, df_dj, df_dk, df_dl, df_dm] = emdlab_flib_gradient5(f, di, dj, dk, dl, dm)
% Computes the gradient of a 5D scalar field f(i,j,k,l,m)
% Inputs:
%   f  - 5D array (n1 × n2 × n3 × n4 × n5)
%   di, dj, dk, dl, dm - grid spacings along each dimension
% Outputs:
%   df_di ... df_dm - partial derivatives along each dimension

[n1, n2, n3, n4, n5] = size(f);

df_di = zeros(size(f));  % ∂f/∂i
df_dj = zeros(size(f));  % ∂f/∂j
df_dk = zeros(size(f));  % ∂f/∂k
df_dl = zeros(size(f));  % ∂f/∂l
df_dm = zeros(size(f));  % ∂f/∂m

for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            for l = 1:n4
                for m = 1:n5

                    % ∂f/∂i
                    if i == 1
                        df_di(i,j,k,l,m) = (f(i+1,j,k,l,m) - f(i,j,k,l,m)) / di;
                    elseif i == n1
                        df_di(i,j,k,l,m) = (f(i,j,k,l,m) - f(i-1,j,k,l,m)) / di;
                    else
                        df_di(i,j,k,l,m) = (f(i+1,j,k,l,m) - f(i-1,j,k,l,m)) / (2*di);
                    end

                    % ∂f/∂j
                    if j == 1
                        df_dj(i,j,k,l,m) = (f(i,j+1,k,l,m) - f(i,j,k,l,m)) / dj;
                    elseif j == n2
                        df_dj(i,j,k,l,m) = (f(i,j,k,l,m) - f(i,j-1,k,l,m)) / dj;
                    else
                        df_dj(i,j,k,l,m) = (f(i,j+1,k,l,m) - f(i,j-1,k,l,m)) / (2*dj);
                    end

                    % ∂f/∂k
                    if k == 1
                        df_dk(i,j,k,l,m) = (f(i,j,k+1,l,m) - f(i,j,k,l,m)) / dk;
                    elseif k == n3
                        df_dk(i,j,k,l,m) = (f(i,j,k,l,m) - f(i,j,k-1,l,m)) / dk;
                    else
                        df_dk(i,j,k,l,m) = (f(i,j,k+1,l,m) - f(i,j,k-1,l,m)) / (2*dk);
                    end

                    % ∂f/∂l
                    if l == 1
                        df_dl(i,j,k,l,m) = (f(i,j,k,l+1,m) - f(i,j,k,l,m)) / dl;
                    elseif l == n4
                        df_dl(i,j,k,l,m) = (f(i,j,k,l,m) - f(i,j,k,l-1,m)) / dl;
                    else
                        df_dl(i,j,k,l,m) = (f(i,j,k,l+1,m) - f(i,j,k,l-1,m)) / (2*dl);
                    end

                    % ∂f/∂m
                    if m == 1
                        df_dm(i,j,k,l,m) = (f(i,j,k,l,m+1) - f(i,j,k,l,m)) / dm;
                    elseif m == n5
                        df_dm(i,j,k,l,m) = (f(i,j,k,l,m) - f(i,j,k,l,m-1)) / dm;
                    else
                        df_dm(i,j,k,l,m) = (f(i,j,k,l,m+1) - f(i,j,k,l,m-1)) / (2*dm);
                    end

                end
            end
        end
    end
end

end
