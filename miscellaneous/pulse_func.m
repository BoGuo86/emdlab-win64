classdef pulse_func < handle
    properties
        t0 = 0;
        D = 0.5;
        min = 0;
        max = 1;
        f = 50;
    end
    methods
        function obj = pulse_func()
        end
        function y = getValue(obj,x)
            T = 1/obj.f;
            x = rem(x-obj.t0,T);
            if x<0
                x = x+T;
            end
            if x < obj.D*T
                y = obj.max;
            else
                y = obj.min;
            end
            
        end
    end
end