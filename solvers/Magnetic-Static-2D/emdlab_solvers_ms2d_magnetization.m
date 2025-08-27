% magnetization object for 2D magnetic-static problems

classdef emdlab_solvers_ms2d_magnetization < handle & matlab.mixin.Copyable

    properties

        % magnetic coercivity [A/m]
        Hc

        % magnetization direction
        magDir

    end

    methods

        function obj = emdlab_solvers_ms2d_magnetization(Hc, magDir)

            if nargin<2
                error('Not enough input arguments.');
            elseif nargin>2
                error('Too many input arguments.');
            end

            if ~isscalar(Hc)
                error('Hc must be an scalar value.')
            end

            obj.Hc = Hc;
            obj.setMagnetizationDirection(magDir);

        end

        function setMagnetizationDirection(obj,value)

            % define direction of magnetization using spatial variables of
            % the global coordinate system
            if isa(value, 'char')

                value = lower(strrep(value,' ',''));
                switch value

                    % uniform magnetization: in the direction of the x-axis
                    case 'x'
                        obj.magDir = [1,0];

                    % uniform magnetization: in the opposite direction of the x-axis
                    case '-x'
                        obj.magDir = [-1,0];

                    % uniform magnetization: in the direction of the y-axis
                    case 'y'
                        obj.magDir = [0,1];

                    % uniform magnetization: in the direction of the y-axis
                    case '-y'
                        obj.magDir = [0,-1];

                    case 'r'
                        obj.magDir = @(x,y) [x,y]/norm([x,y]);

                    case '-r'
                        obj.magDir = @(x,y) -[x,y]/norm([x,y]);

                    case 't'
                        obj.magDir = @(x,y) [y,-x]/norm([x,y]);

                    case '-t'
                        obj.magDir = @(x,y) [-y,x]/norm([x,y]);

                    otherwise
                        error('magDir must be <x>, <y>, <r> or <t>');

                end

            elseif isa(value,'function_handle')

                if nargin(value) ~= 2
                    error('Function handle must get two input <x> and <y> coordinates.')
                end
                obj.magDir = value;

            elseif isa(value,'double')

                [rDim,cDim] = size(value);
                if (rDim~=1)||(cDim~=2)
                    error('Magnetization direction must be a [1x2] row vector.');
                end
                obj.magDir = value;

            else

                error('Magnetization direction must be a <char>, a <function handle> or a <2x1 row vector>.');

            end
        end

        function y = getM(obj,p)

            if isa(obj.magDir,'double')

                y = repmat(-obj.Hc * obj.magDir / norm(obj.magDir), size(p,1), 1);

            elseif isa(obj.magDir,'function_handle')

                y = zeros(size(p,1),2);
                for i = 1:size(p,1)
                    y(i,:) = feval(obj.magDir,p(i,1),p(i,2));
                end

                y = -y * obj.Hc;

            end

        end

        function rotate(obj, varargin)

            if isa(obj.magDir,'double')
                obj.magDir = ext_protate2(obj.magDir, varargin{:});
            end

        end

        function y = getRotate(obj,varargin)
            y = copy(obj);
            if isa(obj.magDir,'double')
                y.magDir = ext_protate2(y.magDir,varargin{:});
            end
        end

    end
end
