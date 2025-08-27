function [f, p] = emdlab_g2d_refinePolylineByRMS(f, p, max_l)

% number of edges
Ne = size(f,1);

for i = 1:Ne
    
    p1 = p(f(i,1),:);
    p2 = p(f(i,2),:);
    
    l = norm(p2-p1);
    
    if l > max_l
        
        n = ceil(l/max_l) + 1;
        x = linspace(p1(1), p2(1), n)';
        y = linspace(p1(2), p2(2), n)';
        
        x = x(2:end-1);
        y = y(2:end-1);
        
        n = n-2;
        
        Np_old = size(p,1);
        
        p = [p;[x,y]];
        
        if n >1
        f_new = [(1:n)+ Np_old;[(2:n)+ Np_old,f(i,2)]]' ;
        else
            f_new = [Np_old+1,f(i,2)];
        end
        
        f(i,2) = Np_old+1;
        f = [f;f_new];
        
    end
    
end



end