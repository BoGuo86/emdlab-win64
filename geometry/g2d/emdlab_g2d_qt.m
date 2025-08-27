classdef emdlab_g2d_qt < handle
    
    properties
        
        coefs (2,3) double;
        
    end
    
    methods
        
        function obj = emdlab_g2d_qt(x1, y1, x2, y2, m1, m2)
            
            A = zeros(6,6);
            A(1,3) = 1;
            A(2,1:3) = 1;
            A(3,6) = 1;
            A(4,4:6) = 1;
            
            if m1 == inf
                A(5,2) = 1;
            else
                A(5,[2,5]) = [m1,-1];
            end
            
            if m2 == inf
                A(6,[1,2]) = [2,1];
            else
                A(6,[1,2,4,5]) = [2*m2,m2,-2,-1];
            end
            
            tmp = A\[x1;x2;y1;y2;0;0];
            
            obj.coefs = tmp([1,2,3;4,5,6]);
            
        end
        
        function plot(obj)
            
            t = linspace(0,1,100);
            
            x = polyval(obj.coefs(1,:),t);
            y = polyval(obj.coefs(2,:),t);
            
            plot(x,y);
        end
        
        function p = getMeshNodes(obj)
            
            t = linspace(0,1,10);
            
            x = polyval(obj.coefs(1,:),t);
            y = polyval(obj.coefs(2,:),t);
            
            p = [x;y]';
            
        end
        
    end
    
end