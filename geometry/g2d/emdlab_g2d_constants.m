classdef emdlab_g2d_constants < handle

    properties (Constant = true, Hidden = true)
        
        % geometry epsilon of leNeh
        gleps = 1e-6;
        % geometry epsilon of angle
        gaeps = 2*atan(1e-6/2);
        
    end
    
end
