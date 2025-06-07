% developer: https://ComProgExpert.com
% 2d magnetic-static coil

classdef emdlab_solvers_ms2d_coil < handle & matlab.mixin.SetGet

    properties
        
        turns (1,1) double {mustBePositive} = 1;
        direction (1,:) char = 'positive';
        parentMatrix (1,:) char = '';
        
    end

    properties (Dependent = true)
        
        sign (1,1) double;
        
    end

    methods

        function obj = emdlab_solvers_ms2d_coil(varargin)
            set(obj, varargin{:});
        end

        function set.turns(obj, value)
            obj.turns = value;
        end

        function set.direction(obj, value)
            value = strrep(value, ' ', '');

            if ismember(lower(value), {'positive', 'negative'})
                obj.direction = value;
            else
                error('Coil direction must be "positive" or "negative".');
            end

        end

        function y = get.sign(obj)

            switch obj.direction
                case 'positive'
                    y = 1;
                case 'negative'
                    y = -1;
            end

        end

    end

end
