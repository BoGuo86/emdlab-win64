classdef emdlab_stp_FACE_BOUND
    properties
        id (1,1) double                 % Unique STEP entity ID
        name (1,:) char = ''            % Optional name string
        edge_loop (1,1) double          % ID of the EDGE_LOOP entity
        orientation (1,1) logical = true % Orientation relative to the face (true = positive sense)
    end

    methods
        function obj = emdlab_stp_FACE_BOUND(id, name, edge_loop, orientation)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 4
                obj.id = id;
                obj.name = name;
                obj.edge_loop = edge_loop;
                obj.orientation = orientation;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', ''''''); % Escape single quotes
            orientationStr = '.T.';
            if ~obj.orientation
                orientationStr = '.F.';
            end
            str = sprintf('#%d=FACE_BOUND(''%s'',#%d,%s);', ...
                obj.id, safeName, obj.edge_loop, orientationStr);
        end
    end
end
