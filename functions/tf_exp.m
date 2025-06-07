classdef tf_exp < handle & matlab.mixin.SetGet
    properties
        amplitude
        power
    end
    methods
        function obj = tf_exp(varargin)
            % default setting
            obj.amplitude = 1;
			obj.power = -1;
			set(obj, varargin{:});
        end
        function y = getValue(obj, time)
            y = obj.amplitude*exp(obj.power*time);
        end
    end
end
