classdef emdlab_stp_CARTESIAN_POINT

    properties
        id (1,1) double; % unique id
        name (1,:) char; % optional name string   
        coordinates (1,3) double % [x,y,z] coordinates 
    end

    methods

        function obj = emdlab_stp_CARTESIAN_POINT(id, name, coordinates)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 3
                obj.id = id;
                obj.name = name;
                obj.coordinates = coordinates;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            % Format coordinates as comma-separated string with high precision
            coordStr = sprintf('%.15g,%.15g,%.15g', obj.coordinates);

            % Escape single quotes in name by doubling them (STEP format)
            safeName = strrep(obj.name, '''', '''''');

            % Build the STEP line
            str = sprintf('#%d = CARTESIAN_POINT(''%s'',(%s));', obj.id, safeName, coordStr);
        end
        
    end
end
