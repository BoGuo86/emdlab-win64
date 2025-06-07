classdef tf_cte < handle & matlab.mixin.SetGet
  properties
    amplitude
  end
  methods
    function obj = tf_cte(varargin)
      % default setting
      obj.amplitude = 1;
      set(obj, varargin{:});
    end
    function y = getValue(obj, time)
      y = obj.amplitude*ones(1, length(time));
    end
  end
end
