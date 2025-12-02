classdef emdlab_stp_SHAPE_REPRESENTATION
    properties
        id (1,1) double                    % Unique STEP entity ID
        name (1,:) char = ''              % Optional name/description
        items (1,:) double = []           % List of IDs of representation items (e.g., solids)
        context (1,1) double              % ID of geometric representation context
    end
    
    methods
        function obj = emdlab_stp_SHAPE_REPRESENTATION(id, name, items, context)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 4
                obj.id = id;
                obj.name = name;
                validateattributes(items, {'double'}, {'row'});
                obj.items = items;
                obj.context = context;
            else
                error('Invalid number of arguments');
            end
        end
        
        function str = toSTEP(obj)
            % Escape single quotes in name
            safeName = strrep(obj.name, '''', '''''');
            % Format items list
            itemsStr = join(string("#" + obj.items), ',');
            % Compose STEP line
            str = sprintf('#%d=SHAPE_REPRESENTATION(''%s'',(%s),#%d);', ...
                obj.id, safeName, itemsStr, obj.context);
        end
    end
end
