function [df_di, df_dj, df_dk, df_dm] = emdlab_flib_gradient4(f, di, dj, dk, dm)

[n1, n2, n3, n4] = size(f);

% Initialize derivative arrays
df_di = zeros(n1, n2, n3, n4);
df_dj = zeros(n1, n2, n3, n4);
df_dk = zeros(n1, n2, n3, n4);
df_dm = zeros(n1, n2, n3, n4);

for ii = 1:n1
    for jj = 1:n2
        for kk = 1:n3
            for mm = 1:n4

                % ∂f/∂i (dim 1)
                if ii == 1
                    df_di(ii,jj,kk,mm) = (f(ii+1,jj,kk,mm) - f(ii,jj,kk,mm)) / di;
                elseif ii == n1
                    df_di(ii,jj,kk,mm) = (f(ii,jj,kk,mm) - f(ii-1,jj,kk,mm)) / di;
                else
                    df_di(ii,jj,kk,mm) = (f(ii+1,jj,kk,mm) - f(ii-1,jj,kk,mm)) / (2*di);
                end

                % ∂f/∂j (dim 2)
                if jj == 1
                    df_dj(ii,jj,kk,mm) = (f(ii,jj+1,kk,mm) - f(ii,jj,kk,mm)) / dj;
                elseif jj == n2
                    df_dj(ii,jj,kk,mm) = (f(ii,jj,kk,mm) - f(ii,jj-1,kk,mm)) / dj;
                else
                    df_dj(ii,jj,kk,mm) = (f(ii,jj+1,kk,mm) - f(ii,jj-1,kk,mm)) / (2*dj);
                end

                % ∂f/∂k (dim 3)
                if kk == 1
                    df_dk(ii,jj,kk,mm) = (f(ii,jj,kk+1,mm) - f(ii,jj,kk,mm)) / dk;
                elseif kk == n3
                    df_dk(ii,jj,kk,mm) = (f(ii,jj,kk,mm) - f(ii,jj,kk-1,mm)) / dk;
                else
                    df_dk(ii,jj,kk,mm) = (f(ii,jj,kk+1,mm) - f(ii,jj,kk-1,mm)) / (2*dk);
                end

                % ∂f/∂m (dim 4)
                if mm == 1
                    df_dm(ii,jj,kk,mm) = (f(ii,jj,kk,mm+1) - f(ii,jj,kk,mm)) / dm;
                elseif mm == n4
                    df_dm(ii,jj,kk,mm) = (f(ii,jj,kk,mm) - f(ii,jj,kk,mm-1)) / dm;
                else
                    df_dm(ii,jj,kk,mm) = (f(ii,jj,kk,mm+1) - f(ii,jj,kk,mm-1)) / (2*dm);
                end

            end
        end
    end
end

end