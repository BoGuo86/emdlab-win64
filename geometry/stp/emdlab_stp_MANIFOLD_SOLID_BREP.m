classdef emdlab_stp_MANIFOLD_SOLID_BREP
    properties
        id (1,1) double                % Unique STEP entity ID
        name (1,:) char = ''           % Optional name
        outer (1,1) double             % ID of CLOSED_SHELL defining the outer shell
    end

    methods
        function obj = emdlab_stp_MANIFOLD_SOLID_BREP(id, name, outer)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 3
                obj.id = id;
                obj.name = name;
                validateattributes(outer, {'double'}, {'scalar'});
                obj.outer = outer;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', '''''');  % Escape single quotes
            str = sprintf('#%d=MANIFOLD_SOLID_BREP(''%s'', #%d);', obj.id, safeName, obj.outer);
        end
    end
end
