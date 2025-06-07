classdef tri_func
    properties
        min = -1;
        max = 1;
        f = 50;
    end
    methods
        function obj = tri_func()
        end
        function y = getValue(obj,x)
            T = 1/obj.f;
            x = rem(x,T);
            if x < T/4
                y = (4*obj.max/T)*x;
            elseif x<3*T/4
                y = (2*(obj.min-obj.max)/T)*(x-T/2);
            else
                y = (4*obj.max/T)*(x-T);
            end
            
        end
    end
end