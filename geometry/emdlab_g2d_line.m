% developer: https://ComProgExpert.com, Ali Jamali-Fard
% 2D segment class

classdef emdlab_g2d_line < handle

    properties

        % a point on the line
        p0 (1,1) emdlab_g2d_point;

        % line direction vector
        u (1,1) emdlab_g2d_point;

    end

    methods
        %% constructor and destructor
        function obj = emdlab_g2d_line(p0, u)

            obj.p0 = p0;
            obj.u = u;
            obj.u.normalize;

        end

        function moveNormal(obj, d)

            n = obj.u.getRotateAroundOrigin(pi/2);
            obj.p0.x = obj.p0.x + d*n.x;
            obj.p0.y = obj.p0.y + d*n.y;

        end

    end

end
