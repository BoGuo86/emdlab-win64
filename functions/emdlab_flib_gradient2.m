% this function calculates the gradients of a f(x,y)

function [df_di, df_dj] = emdlab_flib_gradient2(f, di, dj)

[n1, n2] = size(f);

% Initialize derivative arrays
df_di = zeros(n1, n2);  % ∂f/∂i (along rows)
df_dj = zeros(n1, n2);  % ∂f/∂j (along columns)

for ii = 1:n1
    for jj = 1:n2

        % ∂f/∂i (rows)
        if ii == 1
            df_di(ii,jj) = (f(ii+1,jj) - f(ii,jj)) / di;
        elseif ii == n1
            df_di(ii,jj) = (f(ii,jj) - f(ii-1,jj)) / di;
        else
            df_di(ii,jj) = (f(ii+1,jj) - f(ii-1,jj)) / (2*di);
        end

        % ∂f/∂j (columns)
        if jj == 1
            df_dj(ii,jj) = (f(ii,jj+1) - f(ii,jj)) / dj;
        elseif jj == n2
            df_dj(ii,jj) = (f(ii,jj) - f(ii,jj-1)) / dj;
        else
            df_dj(ii,jj) = (f(ii,jj+1) - f(ii,jj-1)) / (2*dj);
        end

    end
end

end