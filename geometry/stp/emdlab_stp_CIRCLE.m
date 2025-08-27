classdef emdlab_stp_CIRCLE
    properties
        id (1,1) double                         % Unique ID of this entity
        name (1,:) char = ''                   % Optional name
        positionID (1,1) double                % ID of AXIS2_PLACEMENT_3D
        radius (1,1) double {mustBePositive} = 1  % Circle radius
    end

    methods
        function obj = emdlab_stp_CIRCLE(id, name, positionID, radius)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 4
                obj.id = id;
                obj.name = name;
                obj.positionID = positionID;
                obj.radius = radius;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', '''''');  % escape single quotes
            str = sprintf('#%d = CIRCLE(''%s'',#%d,%.15g);', ...
                obj.id, safeName, obj.positionID, obj.radius);
        end
    end
end
