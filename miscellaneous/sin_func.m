classdef sin_func
    properties
        A
        f
        phi
    end
    methods
        function obj = sin_func(A,f,phi)
            obj.A = A;
            obj.f = f;
            obj.phi = phi;
        end
        function y = getValue(obj,x)
            y = obj.A*sin(2*pi*obj.f*x+obj.phi);
        end
    end
end