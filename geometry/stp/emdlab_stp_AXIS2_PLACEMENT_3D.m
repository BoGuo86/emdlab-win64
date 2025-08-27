classdef emdlab_stp_AXIS2_PLACEMENT_3D

    properties
        id (1,1) double;                  % Unique ID
        name (1,:) char = '';             % Optional name
        location (1,1) double;            % ID of CARTESIAN_POINT object
        axis (1,1) double;                % ID of Z axis DIRECTION object
        refDirection (1,1) double;        % ID of X axis DIRECTION object
    end

    methods
        function obj = emdlab_stp_AXIS2_PLACEMENT_3D(id, name, location, axis, refDirection)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 5
                obj.id = id;
                obj.name = name;
                obj.location = location;
                obj.axis = axis;
                obj.refDirection = refDirection;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', '''''');
            str = sprintf('#%d=AXIS2_PLACEMENT_3D(''%s'',#%d,#%d,#%d);', ...
                obj.id, safeName, obj.location, obj.axis, obj.refDirection);
        end
    end
end
