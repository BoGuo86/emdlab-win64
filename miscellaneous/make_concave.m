function x = make_concave(x)

Nx = length(x);

while ~isconcave(x)
    for i = 2:Nx-1
        tmp = (x(i-1) + x(i+1))/2;
        if x(i)>tmp
            x(i) = 1.05*tmp;
        end
    end
end