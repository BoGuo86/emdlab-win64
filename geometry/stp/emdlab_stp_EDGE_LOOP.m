classdef emdlab_stp_EDGE_LOOP
    properties
        id (1,1) double                          % Unique ID for the EDGE_LOOP entity
        name (1,:) char = ''                    % Optional name
        oriented_edges (1,:) double             % Row vector of oriented edge IDs
    end

    methods
        function obj = emdlab_stp_EDGE_LOOP(id, name, oriented_edges)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 3
                obj.id = id;
                obj.name = name;
                validateattributes(oriented_edges, {'double'}, {'row'});
                obj.oriented_edges = oriented_edges;  % Keep as row vector
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', '''''');  % Escape single quotes
            edgeListStr = join(string("#" + obj.oriented_edges), ',');
            str = sprintf('#%d=EDGE_LOOP(''%s'',(%s));', ...
                obj.id, safeName, edgeListStr);
        end
    end
end
