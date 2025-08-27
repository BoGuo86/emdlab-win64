classdef emdlab_stp_ORIENTED_EDGE
    properties
        id (1,1) double                         % Unique STEP entity ID
        name (1,:) char = ''                   % Optional name
        edge_element (1,1) double              % ID of the EDGE_CURVE
        orientation (1,1) logical = true       % Orientation (true = same_sense)
    end

    methods
        function obj = emdlab_stp_ORIENTED_EDGE(id, name, edge_element, orientation)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 4
                obj.id = id;
                obj.name = name;
                obj.edge_element = edge_element;
                obj.orientation = orientation;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', '''''');  % Escape quotes
            orientationStr = '.T.';
            if ~obj.orientation
                orientationStr = '.F.';
            end
            str = sprintf('#%d=ORIENTED_EDGE(''%s'',*,*,#%d,%s);', ...
                obj.id, safeName, obj.edge_element, orientationStr);
        end
    end
end
