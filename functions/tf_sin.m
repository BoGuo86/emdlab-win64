classdef tf_sin < handle & matlab.mixin.SetGet
    properties
        amplitude
        frequency
        phase
    end
    methods
        function obj = tf_sin(varargin)
            % default setting
            obj.amplitude = 1;
			obj.frequency = 50;
			obj.phase = 0;
			set(obj, varargin{:});
        end
        function y = getValue(obj, time)
            y = obj.amplitude*sin(2*pi*obj.frequency*time + obj.phase);
        end
    end
end
