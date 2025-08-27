classdef emdlab_stp_CLOSED_SHELL
    properties
        id (1,1) double                     % Unique STEP entity ID
        name (1,:) char = ''               % Optional name
        faces (1,:) double = []            % List of face IDs forming the closed shell
    end

    methods
        function obj = emdlab_stp_CLOSED_SHELL(id, name, faces)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 3
                obj.id = id;
                obj.name = name;
                validateattributes(faces, {'double'}, {'row'});
                obj.faces = faces;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', '''''');  % Escape any single quotes
            faceStr = join(string("#" + obj.faces), ',');
            str = sprintf('#%d=CLOSED_SHELL(''%s'',(%s));', obj.id, safeName, faceStr);
        end
    end
end
