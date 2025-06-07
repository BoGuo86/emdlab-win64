classdef tf_pulse < handle & matlab.mixin.SetGet
  properties
    ymin (1,1) double = 0;
    ymax (1,1) double = 1;
    period (1,1) double = 0.02;
    duty (1,1) double = 0.5;
    t0 (1,1) double = 0;
  end
  methods
    function obj = tf_pulse(varargin)
      set(obj, varargin{:});
    end
    function y = getValue(obj, t)     
      t = rem(t, obj.period) - obj.t0;
      if t>=0 && t<=obj.duty*obj.period
        y = obj.ymax;
      else
        y = obj.ymin;
      end
    end
  end
end
