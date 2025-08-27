classdef emdlab_stp_ADVANCED_FACE
    properties
        id (1,1) double                     % Unique STEP entity ID
        name (1,:) char = ''               % Optional name
        bounds (1,:) double = []           % Row vector of FACE_BOUND and FACE_OUTER_BOUND IDs
        surface (1,1) double               % ID of the surface (e.g., PLANE, CYLINDRICAL_SURFACE)
        same_sense (1,1) logical = true    % Orientation flag
    end

    methods
        function obj = emdlab_stp_ADVANCED_FACE(id, name, bounds, surface, same_sense)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 5
                obj.id = id;
                obj.name = name;
                validateattributes(bounds, {'double'}, {'row'});
                obj.bounds = bounds;
                obj.surface = surface;
                obj.same_sense = same_sense;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', ''''''); % Escape quotes
            boundsStr = join(string("#" + obj.bounds), ',');
            senseStr = '.T.';
            if ~obj.same_sense
                senseStr = '.F.';
            end
            str = sprintf('#%d=ADVANCED_FACE(%s,(%s),#%d,%s);', ...
                obj.id, safeName,boundsStr, obj.surface, senseStr);
        end
    end
end
