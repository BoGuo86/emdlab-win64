classdef emdlab_stp_VECTOR

    properties
        id (1,1) double; % unique id
        name (1,:) char; % optional name string
        direction (1,1) double % id of DIRECTION object
        magnitude (1,1) double {mustBePositive} = 1;  % magnitude of vector
    end

    methods
        function obj = emdlab_stp_VECTOR(id, name, direction, magnitude)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 4
                obj.id = id;
                obj.name = name;
                obj.direction = direction;
                obj.magnitude = magnitude;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', '''''');
            str = sprintf('#%d=VECTOR(''%s'',#%d,%.15g);', ...
                obj.id, safeName, obj.direction, obj.magnitude);
        end
    end
end
