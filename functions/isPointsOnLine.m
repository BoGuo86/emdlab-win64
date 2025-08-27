function y = isPointsOnLine(p)
    Np = size(p,1);
    if Np<2
        error('Numeber of points at least must be two.');
    end
    if Np == 2
        y = true; return;
    end
    u = p(1,:) - p(2,:);
    for i = 3:Np
        v = p(i,:) - p(1,:);
        if  abs(u(1)*v(2) - u(2)*v(1)) > 1e-6
            y = false; return;
        end
    end
    y = true;
end