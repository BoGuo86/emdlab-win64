classdef emdlab_stp_FACE_OUTER_BOUND
    properties
        id (1,1) double                 % Unique STEP entity ID
        name (1,:) char = ''            % Optional name string
        edge_loop (1,1) double          % ID of the EDGE_LOOP entity
        same_sense (1,1) logical = true % Orientation relative to the face (default true)
    end

    methods
        function obj = emdlab_stp_FACE_OUTER_BOUND(id, name, edge_loop, same_sense)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 3
                obj.id = id;
                obj.name = name;
                obj.edge_loop = edge_loop;
                obj.same_sense = true; % default
            elseif nargin == 4
                obj.id = id;
                obj.name = name;
                obj.edge_loop = edge_loop;
                obj.same_sense = same_sense;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', ''''''); % Escape single quotes
            ss = '.T.';
            if ~obj.same_sense
                ss = '.F.';
            end
            str = sprintf('#%d=FACE_OUTER_BOUND(''%s'',#%d,%s);', ...
                obj.id, safeName, obj.edge_loop, ss);
        end
    end
end
