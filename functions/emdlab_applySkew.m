function y = emdlab_applySkew(x, f, x_skew, N)


y = zeros(1, length(f));

for i = 1:N
    
    y = y + func(x + (i-1)*x_skew/N);
    
end


    function y_tmp = func(x_tmp)
        
        x_tmp = rem(x_tmp, x(end));
        y_tmp = interp1(x, f/N, x_tmp);
        
    end

end