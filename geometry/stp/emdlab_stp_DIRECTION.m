classdef emdlab_stp_DIRECTION

    properties
        id (1,1) double; % unique id
        name (1,:) char; % optional name string
        direction_ratios (1,3) double % direction ratiors representing the unit vector
    end

    methods

        function obj = emdlab_stp_DIRECTION(id, name, direction_ratios)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 3
                obj.id = id;
                obj.name = name;
                normVal = norm(direction_ratios);
                if normVal == 0
                    error('direction_ratios cannot be zero vector.');
                end
                obj.direction_ratios = direction_ratios / normVal;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', '''''');
            dirStr = sprintf('%.15g,%.15g,%.15g', obj.direction_ratios);
            str = sprintf('#%d=DIRECTION(''%s'',(%s));', obj.id, safeName, dirStr);
        end
    end
end
