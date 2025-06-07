function x = make_convex(x)

Nx = length(x);

while ~isconvex(x)
    for i = 2:Nx-1
        tmp = (x(i-1) + x(i+1))/2;
        if x(i)>tmp
            x(i) = 0.95*tmp;
        end
    end
end