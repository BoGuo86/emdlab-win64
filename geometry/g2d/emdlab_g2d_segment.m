% developer: https://ComProgExpert.com, Ali Jamali-Fard
% 2D segment class

classdef emdlab_g2d_segment < handle

    properties

        % initial and desination points
        p0 (1,1) emdlab_g2d_point;
        p1 (1,1) emdlab_g2d_point;

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
        function obj = emdlab_g2d_segment(p0, p1)
            obj.p0 = p0;
            obj.p1 = p1;
        end

        function y = getLength(obj)
            y = getRadialLength(obj.p1 - obj.p0);
        end

        function y = getUnitVector(obj)
            y = getVector(obj.p1 - obj.p0);
            y = y/norm(y);
        end

        %% wireframe mesh functions
        function  nodes = getMeshNodesMinimal(obj)

            nodes = zeros(obj.Nnodes,2);
            nodes(:,1) = [obj.p0.x; obj.p1.x];
            nodes(:,2) = [obj.p0.y; obj.p1.y];

        end

        function nodes = getMeshNodes(obj)

            if obj.isSetNnodes

                nodes = zeros(obj.Nnodes,2);
                nodes(:,1) = linspace(obj.p0.x, obj.p1.x, obj.Nnodes);
                nodes(:,2) = linspace(obj.p0.y, obj.p1.y, obj.Nnodes);

            elseif obj.isSetMaxLength

                Nn = max(ceil(getRadialLength(obj.p1 - obj.p0) / obj.maxLength), 2);
                nodes = zeros(Nn,2);
                nodes(:,1) = linspace(obj.p0.x, obj.p1.x, Nn);
                nodes(:,2) = linspace(obj.p0.y, obj.p1.y, Nn);

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
                u = obj.getUnitVector();
                nodes = zeros(n+1,2);
                nodes(:,1) = [obj.p0.x, obj.p0.x + cumsum(obj.L1 * u(1)* a.^(0:n-2)), obj.p1.x];
                nodes(:,2) = [obj.p0.y, obj.p0.y + cumsum(obj.L1 * u(2)* a.^(0:n-2)), obj.p1.y];

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
            obj.p0.meshSize = mL;
            obj.p1.meshSize = mL;

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

            y = [obj.p0.x + obj.p1.x, obj.p0.y + obj.p1.y]/2;

        end

    end

end