% developer: https://ComProgExpert.com, Ali Jamali-Fard
% 2D point class

classdef emdlab_g3d_point  < handle & matlab.mixin.Copyable & emdlab_g2d_constants

    properties

        % x and y coordinates of the point
        x (1,1) double = 0;
        y (1,1) double = 0;
        z (1,1) double = 0;

        % the object tag
        tag (1,:) char = '';

    end

    methods
        %% constructor: you can define the point with (x,y,z) input or ([x,y,z])
        function obj = emdlab_g3d_point(varargin)

            if nargin == 1

                obj.x = varargin{1}(1);
                obj.y = varargin{1}(2);
                obj.z = varargin{1}(3);

            elseif nargin == 3

                obj.x = varargin{1};
                obj.y = varargin{2};
                obj.z = varargin{3};

            end

        end

        %% set and get coordinates
        function setX(obj, newX)
            obj.x = newX;
        end

        function setY(obj, newY)
            obj.y = newY;
        end

        function setCoordinates(obj, newX, newY)
            obj.x = newX;
            obj.y = newY;
        end

        function val = getX(obj)
            val = obj.x;
        end

        function val = getY(obj)
            val = obj.y;
        end

        % get vector format of the point
        function vec = getVector(obj)

            vec = [obj.x, obj.y];

        end

        % set in vector format ([x,y])
        function setVector(obj, vec)

            obj.x = vec(1);
            obj.y = vec(2);

        end

        function u = getUnitVector(obj)

            u = obj.getVector/obj.getRadialLength;

        end

        %% operation functions
        % add two points
        function p2 = plus(obj, p1)

            p2 = emdlab_g3d_point(obj.x + p1.x, obj.y + p1.y);

        end

        % subtract two point
        function p2 = minus(obj, p1)

            p2 = emdlab_g3d_point(obj.x - p1.x, obj.y - p1.y);

        end

        function p2 = mtimes(ls, rs)

            if isscalar(ls)
                p2 = emdlab_g3d_point(ls*rs.x, ls*rs.y);
            else
                p2 = emdlab_g3d_point(ls.x*rs, ls.y*rs);
            end

        end

        %% distance functions
        % distance of the point from origin
        function y = getRadialLength(obj)
            y = sqrt(obj.x.^2 + obj.y.^2);
        end

        function y = norm(obj)
            y = obj.getRadialLength;
        end

        function y = getDistanceFromOrigin(obj)
            y = sqrt(obj.x.^2 + obj.y.^2);
        end

        %% inner and outer products
        % outer product of two points
        function val = outerProduct(obj, p)
            val = obj.x * p.y - obj.y * p.x;
        end

        function val = cross(obj, p1)
            val = obj.x * p1.y - obj.y * p1.x;
        end

        % inner product of two points
        function val = innerProduct(obj, p)
            val = obj.x * p.x + obj.y * p.y;
        end

        function val = dot(obj, p1)
            val = obj.x * p1.x + obj.y * p1.y;
        end

        %% rotate functions
        function rotateAroundOrigin(obj, angle)

            obj.setVector(ext_protate2(obj.getVector, angle));

        end

        function rotateAroundPoint(obj, p0, angle)

            obj.setVector(ext_protate2(obj.getVector, angle, p0.getVector));

        end

        function p = getRotateAroundOrigin(obj, angle)

            p = emdlab_g3d_point(0,0);
            p.setVector(ext_protate2(obj.getVector, angle));

        end

        function p = getRotateAroundPoint(obj, p0, angle)

            p = emdlab_g3d_point(0,0);
            p.setVector(ext_protate2(obj.getVector, angle, p0.getVector));

        end

        %% shift functions
        function shift(obj, xsh, ysh)

            obj.x = obj.x + xsh;
            obj.y = obj.y + ysh;

        end

        function shiftX(obj, xsh)

            obj.x = obj.x + xsh;

        end

        function shiftY(obj, ysh)

            obj.y = obj.y + ysh;

        end

        function p = getShift(obj, varargin)

            p = copy(obj);
            p.shift(varargin{:});

        end

        function p = getShiftX(obj, varargin)

            p = copy(obj);
            p.shiftX(varargin{:});

        end

        function p = getShiftY(obj, varargin)

            p = copy(obj);
            p.shiftY(varargin{:});

        end

        %% move functions
        function moveX(obj, deltaX)
            obj.x = obj.x + deltaX;
        end

        function moveY(obj, deltaY)
            obj.y = obj.y + deltaY;
        end

        function moveXY(obj, deltaX, deltaY)
            obj.x = obj.x + deltaX;
            obj.y = obj.y + deltaY;
        end

        function moveRadial(obj, deltaR)
            u = obj.getVector / obj.getRadialLength;
            obj.x = obj.x + deltaR * u(1);
            obj.y = obj.y + deltaR * u(2);
        end

        function newObj = getMoveX(obj, varargin)
            newObj = copy(obj);
            newObj.moveX(varargin{:});
        end

        function newObj = getMoveY(obj, varargin)
            newObj = copy(obj);
            newObj.moveY(varargin{:});
        end

        function newObj = getMoveXY(obj, varargin)
            newObj = copy(obj);
            newObj.moveXY(varargin{:});
        end

        function newObj = getMoveRadial(obj, varargin)
            newObj = copy(obj);
            newObj.moveRadial(varargin{:});
        end

    end

end
