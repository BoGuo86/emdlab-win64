classdef emdlab_stp_PLANE < handle
    properties
        id (1,1) double; % unique identifier
        name (1,:) char; % optional name string
        position (1,1) emdlab_stp_AXIS2_PLACEMENT_3D  % reference to AXIS2_PLACEMENT_3D object
    end

    methods
        function obj = emdlab_stp_PLANE(id, name, position)
            if nargin == 0
                return;

            elseif nargin == 1
                obj.id = id;

            elseif nargin == 3
                obj.id = id;
                obj.name = name;
                obj.position = position;
            else
                error('Invalid number of arguments.');
            end
        end

        function str = toSTEP(obj)
            str = sprintf('#%d=PLANE(''%s'',#%d);',  obj.id, obj.name, obj.position.id);
        end
    end
end