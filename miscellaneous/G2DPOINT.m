classdef G2DPOINT
    properties (Constant = true)
        % geometry epsilon for length
        gleps = 1e-6
        % geometry epsilon for angle
        gaeps = 2*atan(1e-6/2)
    end
    properties(SetAccess = private)
        coordinate
    end
    methods
        function obj = G2DPOINT(varargin)
            if nargin<1
                obj.coordinate = [1,0];
            elseif nargin == 1
                obj.coordinate = varargin{1};
            elseif nargin == 2
                obj.coordinate = [varargin{1},varargin{2}];
            end
        end
        function y = getDistance(obj)
            y = norm(obj.coordinate);
        end
        function obj = getRotate(obj,varargin)
            obj.coordinate = ext_protate2(obj.coordinate,varargin{:});
        end
        function obj = getMoveRadial(obj,rsh)
            obj.coordinate = obj.coordinate * (1+rsh/obj.getDistance);
        end
    end
end