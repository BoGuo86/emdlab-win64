classdef tf_cos < handle & matlab.mixin.SetGet
    properties (SetAccess = protected)
        amplitude
        frequency
        phase
    end
    methods
        function obj = tf_cos(varargin)
            % default setting
            obj.amplitude = 1;
			obj.frequency = 50;
			obj.phase = 0;
			set(obj, varargin{:});
        end
        function y = getValue(obj, time)
            y = obj.amplitude*cos(2*pi*obj.frequency*time + obj.phase);
        end
    end
end
