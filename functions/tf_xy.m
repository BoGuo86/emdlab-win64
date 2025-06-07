classdef tf_xy < handle & matlab.mixin.SetGet
    properties
        xyData
        isPeriodic
    end
    methods
        function obj = tf_xy(varargin)
            % default setting
            obj.xyData = [0,0;1,1;2,0];
            obj.isPeriodic = true;
            set(obj, varargin{:});
        end
        function y = getValue(obj, time)
            if obj.isPeriodic
                t_min = min(obj.xyData(:,1));
                t_max = max(obj.xyData(:,1));
                interval = t_max - t_min;
                index = time > t_max;
                time(index) = rem(time(index), interval) + t_min;
				index = time < t_min;
                time(index) = rem(time(index), interval) + t_max;
                y = interp1(obj.xyData(:,1), obj.xyData(:,2), time);
            else
                y = interp1(obj.xyData(:,1), obj.xyData(:,2), time);
                y(isnan(y)) = 0;
            end
        end
        function set.isPeriodic(obj, value)
			if isa(value, 'logical')
				obj.isPeriodic = value;
			else
				error('Input value must be logiacl.');
			end
		end
    end
end
