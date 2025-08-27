classdef emdlab_stp_EDGE_CURVE
    properties
        id (1,1) double;                         % Unique STEP entity ID
        name (1,:) char = '';            % Optional name
        edge_curve (1,1) double;               % ID of curve geometry object (e.g., LINE, CIRCLE)
        start_vertex (1,1) double;              % ID of start vertex (VERTEX_POINT)
        end_vertex (1,1) double;                % ID of end vertex (VERTEX_POINT)
        same_sense (1,1) logical = true;        % TRUE if edge direction matches curve direction
        mesh_size (1,1) double = 2;
    end

    methods
        function obj = emdlab_stp_EDGE_CURVE(id, name, start_vertex, end_vertex, edge_curve, same_sense)
            if nargin == 0
                return;
            elseif nargin == 1
                obj.id = id;
            elseif nargin == 6
                obj.id = id;
                obj.name = name;
                obj.edge_curve = edge_curve;
                obj.start_vertex = start_vertex;
                obj.end_vertex = end_vertex;
                obj.same_sense = same_sense;
            else
                error('Invalid number of arguments.');
            end
        end

        function [x, y, z] = getWireframeMesh(obj, g)
            v1 = g.getEntityHandleByID(obj.start_vertex);
            v2 = g.getEntityHandleByID(obj.end_vertex);
            p1 = g.getEntityHandleByID(v1.vertexPointID);
            p2 = g.getEntityHandleByID(v2.vertexPointID);
            P1 = p1.coordinates;
            P2 = p2.coordinates;

            % Get curve geometry
            try
            curve = g.getEntityHandleByID(obj.edge_curve);
            catch
                x = [];
                y = [];
                z = [];
                return 

            end

            switch class(curve)
                case 'emdlab_stp_LINE'
                    L = norm(P2 - P1);
                    n = max(2, ceil(L / obj.mesh_size));
                    x = linspace(P1(1), P2(1), n);
                    y = linspace(P1(2), P2(2), n);
                    z = linspace(P1(3), P2(3), n);

                case 'emdlab_stp_CIRCLE'
                    % Get AXIS2_PLACEMENT_3D
                    placement = g.getEntityHandleByID(curve.positionID);

                    % Get center point
                    center_point = g.getEntityHandleByID(placement.location);
                    center = center_point.coordinates(:);  % 3x1

                    % Get direction vectors (axis and ref_direction)
                    axis_entity = g.getEntityHandleByID(placement.axis);
                    x_entity = g.getEntityHandleByID(placement.refDirection);
                    z_axis = axis_entity.direction_ratios;         % 3x1 normal to circle plane
                    x_dir  = x_entity.direction_ratios;            % 3x1 reference direction
                    z_axis = z_axis / norm(z_axis);            % Ensure unit vectors
                    x_dir = x_dir / norm(x_dir);
                    y_dir = cross(z_axis, x_dir);              % Construct local y-dir

                    % Vectors from center to start and end points
                    v1 = P1(:) - center;
                    v2 = P2(:) - center;

                    % Project to local 2D plane
                    u1 = [dot(v1, x_dir); dot(v1, y_dir)];
                    u2 = [dot(v2, x_dir); dot(v2, y_dir)];

                    % Local angles
                    theta1 = atan2(u1(2), u1(1));
                    theta2 = atan2(u2(2), u2(1));

                    % Adjust angle direction based on same_sense
                    if obj.same_sense
                        if theta2 <= theta1
                            theta2 = theta2 + 2*pi;
                        end
                    else
                        if theta1 <= theta2
                            theta1 = theta1 + 2*pi;
                        end
                    end

                    % Sampling arc
                    angle_span = abs(theta2 - theta1);
                    arc_length = curve.radius * angle_span;
                    n = max(2, ceil(arc_length / obj.mesh_size));
                    theta = linspace(theta1, theta2, n);
                    xy_local = curve.radius * [cos(theta); sin(theta)];

                    % Map local 2D points back to 3D
                    points3D = center + x_dir' .* repmat(xy_local(1,:),3,1) + y_dir' .* repmat(xy_local(2,:),3,1);
                    x = points3D(1,:);
                    y = points3D(2,:);
                    z = points3D(3,:);
                otherwise
                    error('Unsupported curve type: %s', class(curve));
            end

        end

        function str = toSTEP(obj)
            safeName = strrep(obj.name, '''', '''''');  % Escape single quotes
            senseStr = '.T.';
            if obj.same_sense == false
                senseStr = '.F.';
            end
            str = sprintf('#%d=EDGE_CURVE(''%s'',#%d,#%d,#%d,%s);', ...
                obj.id, safeName, obj.start_vertex, obj.end_vertex, obj.edge_curve, senseStr);
        end
    end
end
