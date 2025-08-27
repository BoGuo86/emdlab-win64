% two-dimensional magnetic-static coil
% coil is made of a series of coil arms: there is no parallel branches

classdef emdlab_solvers_ms2d_coil < handle & matlab.mixin.SetGet

    properties

        % the number of coil turns
        coilArms (1,:) string;

        % coil index
        ci (1,1) double;

        % coil arm directions: '+1' -> z+, or '-1' -> z-
        directions (1,:) double;

        % current
        current (1,:) double = 0;

        % flux linkage
        fluxLinkage (1,:) double = 0;

        % dc resistance of the coil
        Rdc (1,1) double;      

        % length of the end winding
        Lend (1,1) double {mustBeNonnegative} = 0;

    end

    properties (Dependent = true)

        % number of coil arms
        NcoilArms (1,1) double;

    end

    methods

        function obj = emdlab_solvers_ms2d_coil()
        end   

        function addCoilArm(obj, newCoilArmName, direction)
            obj.coilArms(end+1) = newCoilArmName;
            obj.directions(end+1) = direction;
        end

        function y = get.NcoilArms(obj)
            y = numel(obj.coilArms);
        end

        function setCurrent(obj, value)
            if (isscalar(value) && isnumeric(value))
                obj.current = value;
            else
                error('Coil current must be a scalar numeric value.');
            end
        end

    end
end
