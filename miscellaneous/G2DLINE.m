classdef G2DLINE
    properties
        p1
        p2
    end
    methods
        function obj = G2DLINE(p1,p2)
            obj.p1 = p1;
            obj.p2 = p2;
        end
        function y = getu(obj)
            y = obj.p2 - obj.p1;
            y = y/norm(y);
        end
        function y = getn(obj)
            y = obj.getu;
            y = [y(2),-y(1)];
        end
        function y = getd(obj,p)
            y = obj.getn;
            y = p*y' - obj.p1*y';
        end
    end
    
end