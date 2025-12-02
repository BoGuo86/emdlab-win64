% EMDLAB: Electrical Machines Design Laboratory
% 2D segment class

classdef emdlab_g2d_ray < handle

    properties

        % starting point of the ray
        p0 (1,1) emdlab_g2d_point = emdlab_g2d_point(0,0);

        % unit vector in the direction of the ray
        u (1,1) emdlab_g2d_point = emdlab_g2d_point(1,0);

    end

    methods
        %% constructor and destructor
        function obj = emdlab_g2d_ray(varargin)

            if nargin
                obj.p0 = varargin{1};
                obj.u = varargin{2};
                obj.u.normalize;
            end

        end

        function moveNormal(obj, d)

            n = obj.u.getRotateAroundOrigin(pi/2);
            obj.p0.x = obj.p0.x + d*n.x;
            obj.p0.y = obj.p0.y + d*n.y;

        end

        function p = getPoint(obj, t)
            p = obj.p0 + t * obj.u;
        end

        function location = classifyPoint(obj, p)
            arguments
                obj
                p (1,1) emdlab_g2d_point
            end
            tmp = obj.u.cross(p - obj.p0);
            if tmp>0
                location = emdlab_g2d_location.pl_left;
            elseif tmp<0
                location = emdlab_g2d_location.pl_right;
            else
                location = emdlab_g2d_location.pl_on;
            end

        end

    end

end
