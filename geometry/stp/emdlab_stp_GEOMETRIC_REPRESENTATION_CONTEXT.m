classdef emdlab_stp_GEOMETRIC_REPRESENTATION_CONTEXT
    properties
        id (1,1) double               % Unique STEP entity ID
        name (1,:) char = ''          % Optional name or description
        dimension (1,1) double        % Dimension (e.g., 2 or 3)
        coordinateSpaceDimension (1,1) double = NaN % Optional coordinate space dimension
    end

    methods
        function obj = emdlab_stp_GEOMETRIC_REPRESENTATION_CONTEXT(id, name, dimension, coordinateSpaceDimension)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 3
                obj.id = id;
                obj.name = name;
                obj.dimension = dimension;
                obj.coordinateSpaceDimension = dimension; % default to same
            elseif nargin == 4
                obj.id = id;
                obj.name = name;
                obj.dimension = dimension;
                obj.coordinateSpaceDimension = coordinateSpaceDimension;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', '''''');
            if isnan(obj.coordinateSpaceDimension)
                % Omit coordinateSpaceDimension if not provided
                str = sprintf('#%d=GEOMETRIC_REPRESENTATION_CONTEXT(''%s'', %d);', ...
                    obj.id, safeName, obj.dimension);
            else
                str = sprintf('#%d=GEOMETRIC_REPRESENTATION_CONTEXT(''%s'', %d, %d);', ...
                    obj.id, safeName, obj.dimension, obj.coordinateSpaceDimension);
            end
        end
    end
end
