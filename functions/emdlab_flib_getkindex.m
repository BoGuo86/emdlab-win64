function y = emdlab_flib_getkindex(dim)

y = zeros(dim,dim);

index = 0;
for i = 1:dim
    for j = 1:i
        index = index + 1;
        y(i,j) = index;
    end
end

y = y + tril(y,-1)';
y = y(:);

end