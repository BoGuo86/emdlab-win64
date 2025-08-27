classdef emdlab_stp_VERTEX
    properties
        id (1,1) double; % unique id
        name (1,:) char = ''; % optional name string 
        vertexPointID (1,1) double; % ID of associated CARTESIAN_POINT
    end

    methods

        function obj = emdlab_stp_VERTEX(id, name, vertexPointID)
            if nargin == 0
                return;
            elseif nargin == 3
                obj.id = id;
                obj.name = name;
                obj.vertexPointID = vertexPointID;
            elseif nargin == 1
                obj.id = id;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', '''''');  % Escape single quotes
            str = sprintf('#%d = VERTEX_POINT(''%s'',#%d);', ...
                obj.id, safeName, obj.vertexPointID);
        end
        
    end
end
