% EMDLAB: Electrical Machines Design Laboratory
% 2d spline class

classdef emdlab_g2d_spline < handle

    properties

        % spline points
        pts (1,:) emdlab_g2d_point;

        % splines for x and y coordinates
        sx (1,1);
        sy (1,1);

        % properties needed for wireframe mesh
        Nnodes (1,1) double {mustBeInteger, mustBeNonnegative} = 2;
        isSetNnodes (1,1) logical = true;
        maxLength (1,1) double {mustBeNonnegative} = 1;
        isSetMaxLength (1,1) logical = false;
        L1 (1,1) double {mustBeNonnegative} = 1;
        L2 (1,1) double {mustBeNonnegative} = 1;
        isSetL1L2 (1,1) logical = false;

        % object tag
        tag (1,:) char;

    end

    methods
        %% constructor and destructor
        function obj = emdlab_g2d_spline(pts)

            obj.pts = pts;
            Np = numel(pts);
            x = zeros(1,Np);
            y = zeros(1,Np);

            for i = 1:Np
                x(i) = pts(i).x;
                y(i) = pts(i).y;
            end

            t = [0, sqrt(diff(x).^2 + diff(y).^2)];
            t = cumsum(t) / sum(t);

            obj.sx = spline(t, x);
            obj.sy = spline(t, y);

        end

        function y = getLength(obj)
            xi = ppval(obj.sx, linspace(0, 1, 1000));
            yi = ppval(obj.sy, linspace(0, 1, 1000));
            y = [diff(xi);diff(yi)];
            y = sum(sqrt(sum(y.^2)));
        end

        function y = getUnitVector(obj)
            y = getVector(obj.p1 - obj.pts);
            y = y/norm(y);
        end

        %% wireframe mesh functions
        function  nodes = getMeshNodesMinimal(obj)

            nodes = zeros(20,2);
            nodes(:,1) = ppval(obj.sx, linspace(0, 1, 20));
            nodes(:,2) = ppval(obj.sy, linspace(0, 1, 20));

        end

        function nodes = getMeshNodes(obj)

            if obj.isSetNnodes

                nodes = zeros(obj.Nnodes,2);
                nodes(:,1) = ppval(obj.sx, linspace(0, 1, obj.Nnodes));
                nodes(:,2) = ppval(obj.sy, linspace(0, 1, obj.Nnodes));

            elseif obj.isSetMaxLength

                Nn = max(ceil(obj.getLength / obj.maxLength), 2);
                nodes = zeros(Nn,2);
                nodes(:,1) = ppval(obj.sx, linspace(0, 1, Nn));
                nodes(:,2) = ppval(obj.sy, linspace(0, 1, Nn));

            else

                % suppose we have n cords, n+1 points
                % [L1, a*L1, a^2*L1, ...., a^(n-1)*L1], a^(n-1)*L1 = L2
                L = obj.getLength;
                n = 2;
                a = nthroot(obj.L2/obj.L1,n-1);
                errOld = abs(L/obj.L1 - sum(a.^(0:n-1)));
                while true
                    n = n + 1;
                    a = nthroot(obj.L2/obj.L1,n-1);
                    errNew = abs(L/obj.L1 - sum(a.^(0:n-1)));
                    if errNew > errOld
                        n = n - 1;
                        break;
                    end
                    errOld = errNew;
                end

                err_fcn = @(x) sum(x.^(0:n-1))-L/obj.L1;
                a = fzero(err_fcn,1);
                nodes = zeros(n+1,2);
                a = cumsum(obj.L1 * a.^(0:n-2)/L);
                nodes(:,1) = [obj.pts(1).x, ppval(obj.sx, a), obj.pts(end).x];
                nodes(:,2) = [obj.pts(1).y, ppval(obj.sy, a), obj.pts(end).y];                

            end

        end

        function setNnodes(obj, n)

            obj.isSetMaxLength = false;
            obj.isSetNnodes = true;
            obj.isSetL1L2 = false;

            obj.Nnodes = n;

        end

        function setMaxLength(obj, mL)

            obj.isSetMaxLength = true;
            obj.isSetNnodes = false;
            obj.isSetL1L2 = false;

            obj.maxLength = mL;

        end

        function setL1L2(obj, newL1, newL2)

            if (newL1 + newL2) > obj.getLength()
                obj.setNnodes(2);
                return;
            end

            obj.isSetMaxLength = false;
            obj.isSetNnodes = false;
            obj.isSetL1L2 = true;

            obj.L1 = newL1;
            obj.L2 = newL2;

        end

        function y = getCenter(obj)

            y = [ppval(obj.sx, 0.5), ppval(obj.sy, 0.5)];

        end

    end

end