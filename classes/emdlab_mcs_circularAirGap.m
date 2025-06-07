% developer: https://ComProgExpert.com, Ali Jamali-Fard
% moving contacts -> circular air-gap
% it used for full air-gap modeling of rotating electrical machines
% this object is for meshing of a circular band

classdef emdlab_mcs_circularAirGap < handle & emdlab_g2d_constants & matlab.mixin.SetGet

    properties

        % mesh: triangular mesh zone
        m;

        % inner points
        ips (:,2) double;

        % outer points
        ops (:,2) double;

        % inner point angles
        ipas (:,1) double;

        % outer point angles
        opas (:,1) double;

        % minimum mesh length
        h0 (1,1) double;

        % number of layers for mesh
        Nlayer (1,1) double {mustBeInteger, mustBePositive} = 1;

        % inner circle center
        icc (1,2) double = [0, 0];

        % outer circle center
        occ (1,2) double = [0, 0];

        % number of inner points
        Nips (1,1) double;

        % number of outer points
        Nops (1,1) double;

        % inner circle radius
        rin (1,1) double;

        % outer circle radius
        rout (1,1) double;

    end

    methods

        function obj = emdlab_mcs_circularAirGap(ips, ops, Nlayer, varargin)

            obj.ips = ips;
            obj.ops = ops;
            
            obj.Nips = size(ips, 1);
            obj.Nops = size(ops, 1);

            obj.rin = obj.get_rin;
            obj.rout = obj.get_rout;
            
            obj.Nlayer = Nlayer;

            % icc and occ must be set manually if they are not [0,0]
%             set(obj, varargin{:});
            % chekers: inner points must be insode outer points
            obj.sortPoints;
            % evaluation of h0
            l1 = sqrt(sum((obj.ips - circshift(obj.ips, -1)).^2, 2));
            l2 = sqrt(sum((obj.ops - circshift(obj.ops, -1)).^2, 2));
            obj.h0 = mean([l1; l2]);
            obj.updateMesh;
        end

        function y = get_rin(obj)
            % calculation of inner circle radius
            tmp = [obj.ips(:, 1) - obj.icc(1), obj.ips(:, 2) - obj.icc(2)];
            tmp = sqrt(sum(tmp.^2, 2));
            y = mean(tmp);

            if sum(abs(tmp - y)) > obj.Nips * obj.gleps
                error('Inner points do not form a cirlce.');
            end

        end

        function y = get_rout(obj)
            % calculation of inner circle radius
            tmp = [obj.ops(:, 1) - obj.occ(1), obj.ops(:, 2) - obj.occ(2)];
            tmp = sqrt(sum(tmp.^2, 2));
            y = mean(tmp);

            if sum(abs(tmp - y)) > obj.Nops * obj.gleps
                error('Outer points do not form a cirlce.');
            end

        end

        function sortPoints(obj)
            % inner points
            obj.ipas = atan_02pi([obj.ips(:, 1) - obj.icc(1), obj.ips(:, 2) - obj.icc(2)]);
            [~, index] = sort(obj.ipas);
            obj.ips = obj.ips(index, :);
            obj.ipas = obj.ipas(index, :);
            % outer points
            obj.opas = atan_02pi([obj.ops(:, 1) - obj.occ(1), obj.ops(:, 2) - obj.occ(2)]);
            [~, index] = sort(obj.opas);
            obj.ops = obj.ops(index, :);
            obj.opas = obj.opas(index, :);
        end

        function updateMesh(obj)
            obj.checkInside;
            ces1 = [1:obj.Nips; [2:obj.Nips, 1]]';
            ces2 = [1:obj.Nops; [2:obj.Nops, 1]]';

            if obj.Nlayer > 1
                rmid = linspace(obj.rin, obj.rout, obj.Nlayer + 1);
                rmid(1) = [];
                rmid(end) = [];
                Nr = length(rmid);
                n = ceil(pi * (obj.rin + obj.rout) / obj.h0);
                teta = linspace(0, 2 * pi, n + 1)';
                teta(end) = [];
                xp = cos(teta) * rmid;
                yp = sin(teta) * rmid;
                e1 = reshape((1:n * Nr)', [], Nr);
                e2 = circshift(e1, -1);
                ces3 = [e1(:), e2(:)] + obj.Nips + obj.Nops;
            else
                xp = [];
                yp = [];
                ces3 = zeros([], 2);
            end

            ces = [fliplr(ces1); ces2 + obj.Nips];
            obj.m = emdlab_m2d_mm(ces, [obj.ips; obj.ops; xp(:), yp(:)], ces3);
        end

        function checkInside(obj)
            tmp = [obj.ips(:, 1) - obj.occ(1), obj.ips(:, 1) - obj.occ(1)];
            tmp = sqrt(sum(tmp.^2, 2));

            if any(abs(tmp - obj.rout) < obj.gleps)
                error('Inner cirlce must be completely inside ouetr circle.');
            end

        end

        function obj = rotateInnerStatic(obj, alpha)
            obj.ips = ext_protate2(obj.ips, alpha, obj.icc);
            obj.updateMesh;
        end
        
        function obj = rotateOuterStatic(obj, alpha)
            obj.ops = ext_protate2(obj.ops, alpha, obj.occ);
            obj.updateMesh;
        end
        
        function obj = rotateInner(obj, varargin)
            obj.ips = ext_protate2(obj.ips, varargin{:});
            obj.updateMesh;
        end

        function obj = shiftInner(obj, varargin)
            obj.ips = ext_pshift2(obj.ips, varargin{:});
            obj.icc = ext_pshift2(obj.icc, varargin{:});
            obj.updateMesh;
        end

        function obj = rotateOuter(obj, varargin)
            obj.ops = ext_protate2(obj.ops, varargin{:});
            obj.updateMesh;
        end

        function obj = shiftOuter(obj, varargin)
            obj.ops = ext_pshift2(obj.ops, varargin{:});
            obj.occ = ext_pshift2(obj.occ, varargin{:});
            obj.updateMesh;
        end

    end

end
