function z = emdlab_g2d_3parc(x1, y1, x2, y2, x3, y3, varargin)

y = fsolve(@sub_func, [1,1,1]);
z = emdlab_g2d_arc(emdlab_g2d_point(y(1), y(2)), emdlab_g2d_point(x1,y1), emdlab_g2d_point(x3,y3), varargin{:});

    function e = sub_func(x)
        
        e(1) = (x(1) - x1)^2 + (x(2) - y1)^2 - x(3)^2;
        e(2) = (x(1) - x2)^2 + (x(2) - y2)^2 - x(3)^2;
        e(3) = (x(1) - x3)^2 + (x(2) - y3)^2 - x(3)^2;
        
    end

end