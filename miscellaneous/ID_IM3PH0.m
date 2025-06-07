classdef ID_IM3PH0
    properties
        Rso = 50;
        D = 0.6;
        ns = 36;
        nr = 24;
    end
    properties (Dependent=true)
        Rro
        Tr
        Ts
    end
    methods
        function obj = ID_IM3PH0()
        end
        function y = get.Rro(obj)
            y = obj.Rso*obj.D;
        end
        function y = get.Ts(obj)
            y = 2*pi/obj.ns;
        end
        function y = get.Tr(obj)
            y = 2*pi/obj.nr;
        end
    end
    
end