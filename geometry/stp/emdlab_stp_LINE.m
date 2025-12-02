classdef emdlab_stp_LINE
    properties
        id (1,1) double                         % Unique STEP ID
        name (1,:) char = ''                   % Optional name
        position (1,1) double                  % ID of AXIS2_PLACEMENT_3D
        direction (1,1) double                 % ID of DIRECTION
    end

    methods
        function obj = emdlab_stp_LINE(id, name, position, direction)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 4
                obj.id = id;
                obj.name = name;
                obj.position = position;
                obj.direction = direction;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', '''''');  % Escape single quotes
            str = sprintf('#%d=LINE(''%s'',#%d,#%d);', ...
                obj.id, safeName, obj.position, obj.direction);
        end
    end
end
