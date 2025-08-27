function [df_di, df_dj, df_dk, df_dl, df_dm, df_dn] = emdlab_flib_gradient6(f, di, dj, dk, dl, dm, dn)
% Computes the numerical gradients of a 6D scalar field f(i,j,k,l,m,n)
% Inputs:
%   f   - 6D array
%   di, dj, dk, dl, dm, dn - spacing along each dimension
% Outputs:
%   df_di ... df_dn - partial derivatives along each axis

[n1, n2, n3, n4, n5, n6] = size(f);

df_di = zeros(size(f));
df_dj = zeros(size(f));
df_dk = zeros(size(f));
df_dl = zeros(size(f));
df_dm = zeros(size(f));
df_dn = zeros(size(f));

for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            for l = 1:n4
                for m = 1:n5
                    for n = 1:n6

                        % ∂f/∂i
                        if i == 1
                            df_di(i,j,k,l,m,n) = (f(i+1,j,k,l,m,n) - f(i,j,k,l,m,n)) / di;
                        elseif i == n1
                            df_di(i,j,k,l,m,n) = (f(i,j,k,l,m,n) - f(i-1,j,k,l,m,n)) / di;
                        else
                            df_di(i,j,k,l,m,n) = (f(i+1,j,k,l,m,n) - f(i-1,j,k,l,m,n)) / (2*di);
                        end

                        % ∂f/∂j
                        if j == 1
                            df_dj(i,j,k,l,m,n) = (f(i,j+1,k,l,m,n) - f(i,j,k,l,m,n)) / dj;
                        elseif j == n2
                            df_dj(i,j,k,l,m,n) = (f(i,j,k,l,m,n) - f(i,j-1,k,l,m,n)) / dj;
                        else
                            df_dj(i,j,k,l,m,n) = (f(i,j+1,k,l,m,n) - f(i,j-1,k,l,m,n)) / (2*dj);
                        end

                        % ∂f/∂k
                        if k == 1
                            df_dk(i,j,k,l,m,n) = (f(i,j,k+1,l,m,n) - f(i,j,k,l,m,n)) / dk;
                        elseif k == n3
                            df_dk(i,j,k,l,m,n) = (f(i,j,k,l,m,n) - f(i,j,k-1,l,m,n)) / dk;
                        else
                            df_dk(i,j,k,l,m,n) = (f(i,j,k+1,l,m,n) - f(i,j,k-1,l,m,n)) / (2*dk);
                        end

                        % ∂f/∂l
                        if l == 1
                            df_dl(i,j,k,l,m,n) = (f(i,j,k,l+1,m,n) - f(i,j,k,l,m,n)) / dl;
                        elseif l == n4
                            df_dl(i,j,k,l,m,n) = (f(i,j,k,l,m,n) - f(i,j,k,l-1,m,n)) / dl;
                        else
                            df_dl(i,j,k,l,m,n) = (f(i,j,k,l+1,m,n) - f(i,j,k,l-1,m,n)) / (2*dl);
                        end

                        % ∂f/∂m
                        if m == 1
                            df_dm(i,j,k,l,m,n) = (f(i,j,k,l,m+1,n) - f(i,j,k,l,m,n)) / dm;
                        elseif m == n5
                            df_dm(i,j,k,l,m,n) = (f(i,j,k,l,m,n) - f(i,j,k,l,m-1,n)) / dm;
                        else
                            df_dm(i,j,k,l,m,n) = (f(i,j,k,l,m+1,n) - f(i,j,k,l,m-1,n)) / (2*dm);
                        end

                        % ∂f/∂n
                        if n == 1
                            df_dn(i,j,k,l,m,n) = (f(i,j,k,l,m,n+1) - f(i,j,k,l,m,n)) / dn;
                        elseif n == n6
                            df_dn(i,j,k,l,m,n) = (f(i,j,k,l,m,n) - f(i,j,k,l,m,n-1)) / dn;
                        else
                            df_dn(i,j,k,l,m,n) = (f(i,j,k,l,m,n+1) - f(i,j,k,l,m,n-1)) / (2*dn);
                        end

                    end
                end
            end
        end
    end
end
end
