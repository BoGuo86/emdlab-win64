classdef emdlab_stp_MODEL < handle

    properties
        entities = {};   % Cell array to store entities
        ids (1,:) double; % available ids
    end

    methods

        function id = getAUniqueID(obj)
            if isempty(obj.ids)
                id = 1;
            else
                id = max(obj.ids) + 1;
            end
        end

        function entityHandle = getEntityHandleByID(obj, entityID)
            for i = 1:numel(obj.entities)
                e = obj.entities{i};
                if e.id == entityID
                    entityHandle = e;
                    return
                end
            end
            error('Entity was not found.');
        end

        function id = checkIDExistance(obj, id)
            if ~ismember(id, obj.ids)
                error('The id is not assignd to another entity.');
            end
        end

        function id = checkIDNonExistance(obj, id)
            if ismember(id, obj.ids)
                error('The id is already assignd to another entity.');
            end
        end

        function [id, entity] = addEntityByDefaultID(obj, entityName, varargin)
            id = obj.getAUniqueID;
            entity = eval(sprintf('emdlab_stp_%s(id, varargin{:});', entityName));
            obj.ids(end+1) = id;
            obj.entities{end+1} = entity;
        end

        function [id, entity] = addEntity(obj, entityName, varargin)
            id = obj.checkIDNonExistance(varargin{1});
            entity = eval(sprintf('emdlab_stp_%s(varargin{:});', entityName));
            obj.ids(end+1) = id;
            obj.entities{end+1} = entity;
        end

        function [id, entity] = addEntityByText(obj, str)
            entityName = extractBetween(str,'=','(');
            entityName = char(erase(entityName, ' '));
            str = strsplit(str,'=');
            id = str2double(str{1}(2:end));
            str = string(str{2});
            str = erase(str, ' ');
            switch entityName
                case 'CARTESIAN_POINT'
                    str = extractBetween(str,',(','))');
                    str = strsplit(str,',');
                    if numel(str) == 2, str{3} = '0'; end
                    [~,entity] = obj.addEntity('CARTESIAN_POINT',id,'', [str2double(str{1}),str2double(str{2}),str2double(str{3})]);
                case 'DIRECTION'
                    str = extractBetween(str,',(','))');
                    str = strsplit(str,',');
                    if numel(str) == 2, str{3} = '0'; end
                    [~,entity] = obj.addEntity('DIRECTION',id,'', [str2double(str{1}),str2double(str{2}),str2double(str{3})]);
                case 'VERTEX'
                    str = extractBetween(str,'(',')');
                    str = strsplit(str,',');
                    index = str2double(str{2}(2:end));
                    [~,entity] = obj.addEntity('VERTEX',id,'',index);
                case 'VERTEX_POINT'
                    str = extractBetween(str,'(',')');
                    str = strsplit(str,',');
                    index = str2double(str{2}(2:end));
                    [~,entity] = obj.addEntity('VERTEX',id,'',index);
                case 'AXIS2_PLACEMENT_3D'
                    str = extractBetween(str,'(',')');
                    str = strsplit(str,',');
                    [~,entity] = obj.addEntity('AXIS2_PLACEMENT_3D',id,'',str2id(str{2}),str2id(str{3}),str2id(str{4}));
                case 'LINE'
                    str = extractBetween(str,'(',')');
                    str = strsplit(str,',');
                    [~,entity] = obj.addEntity('LINE',id,'',str2id(str{2}),str2id(str{3}));
                case 'CIRCLE'
                    str = extractBetween(str,'(',')');
                    str = strsplit(str,',');
                    [~,entity] = obj.addEntity('CIRCLE',id,'',str2id(str{2}),str2double(str{3}));
                case 'EDGE_CURVE'
                    str = extractBetween(str,'(',')');
                    str = strsplit(str,',');
                    [~,entity] = obj.addEntity('EDGE_CURVE',id,'',str2id(str{2}),str2id(str{3}),str2id(str{4}),str2bool(str{5}));
                case 'ORIENTED_EDGE'
                    str = extractBetween(str,'(',')');
                    str = strsplit(str,',');
                    [~,entity] = obj.addEntity('ORIENTED_EDGE',id,'',str2id(str{4}),str2bool(str{5}));
                case 'EDGE_LOOP'
                    str = extractBetween(str,',(','))');
                    str = strsplit(str,',');
                    tmp = zeros(1,numel(str));
                    for i = 1:numel(str)
                        tmp(i) = str2id(str{i});
                    end
                    [~,entity] = obj.addEntity('EDGE_LOOP',id,'',tmp);
                case 'FACE_OUTER_BOUND'
                    str = extractBetween(str,'(',')');
                    str = strsplit(str,',');
                    [~,entity] = obj.addEntity('FACE_OUTER_BOUND',id,'',str2id(str{2}),str2bool(str{3}));
                case 'FACE_BOUND'
                    str = extractBetween(str,'(',')');
                    str = strsplit(str,',');
                    [~,entity] = obj.addEntity('FACE_BOUND',id,'',str2id(str{2}),str2bool(str{3}));
                case 'ADVANCED_FACE'
                    strBounds = extractBetween(str,',(','),');
                    strBounds = strsplit(strBounds,',');
                    tmp = zeros(1,numel(strBounds));
                    for i = 1:numel(strBounds)
                        tmp(i) = str2id(strBounds{i});
                    end
                    str = regexprep(str, '\(\#\d+(,\#\d+)*\)', '*');
                    str = strsplit(str,',');
                    [~,entity] = obj.addEntity('ADVANCED_FACE ',id,'',tmp,str2id(str{3}),str2bool(str{4}));
            end

            function y = str2id(str)
                y = str2double(str(2:end));
            end

            function y = str2bool(str)
                if strcmpi(str,'.T.')
                    y = true;
                else
                    y = false;
                end
            end

        end


        function showModel(obj)
            [f,ax] = emdlab_r3d_mesh;
            hold on;
            p = zeros([],3);
            nodes = {};
            cl = {};
            tmpIndex = 0;
            for i = 1:numel(obj.entities)
                e = obj.entities{i};
                if isa(e, 'emdlab_stp_VERTEX')
                    point = obj.getEntityHandleByID(e.vertexPointID);
                    p(end+1,:) = point.coordinates;
                end
                if isa(e, 'emdlab_stp_EDGE_CURVE')
                    [x,y,z] = e.getWireframeMesh(obj);
                    nodes{end+1} = [x;y;z];
                    Np = size(nodes{end},2);
                    cl{end+1} = [1:(Np-1);2:Np] + tmpIndex;
                    tmpIndex = tmpIndex + Np;
                end
            end
            plot3(p(:,1),p(:,2),p(:,3),'s', 'parent', ax, 'MarkerEdgeColor','k', 'LineWidth',2);

            nodes = cell2mat(nodes)';
            cl = cell2mat(cl)';
            patch('faces', cl, 'vertices', nodes, 'edgecolor', 'b', 'linewidth',1.2);

            set(ax,'clipping','off')
            %             zoom(ax,'on');
            axis off equal;
            set(f,'Visible','on');
        end

        function addRectangle(obj, x0, y0, w, h)


            o = obj.addEntityByDefaultID('CARTESIAN_POINT','',[0,0,0]);
            xp = obj.addEntityByDefaultID('DIRECTION','',[1,0,0]);
            yp = obj.addEntityByDefaultID('DIRECTION','',[0,1,0]);
            xn = obj.addEntityByDefaultID('DIRECTION','',[-1,0,0]);
            yn = obj.addEntityByDefaultID('DIRECTION','',[0,-1,0]);
            z = obj.addEntityByDefaultID('DIRECTION','',[0,0,1]);
            c = obj.addEntityByDefaultID('AXIS2_PLACEMENT_3D','',o,z,xp);
            xy_plane = obj.addEntityByDefaultID('PLANE','',c);

            p1 = obj.addEntityByDefaultID('CARTESIAN_POINT','',[x0,y0,0]);
            p2 = obj.addEntityByDefaultID('CARTESIAN_POINT','',[x0+w,y0,0]);
            p3 = obj.addEntityByDefaultID('CARTESIAN_POINT','',[x0+w,y0+h,0]);
            p4 = obj.addEntityByDefaultID('CARTESIAN_POINT','',[x0,y0+h,0]);

            v1 = obj.addEntityByDefaultID('VECTOR','',xp,1);
            v2 = obj.addEntityByDefaultID('VECTOR','',yp,1);
            v3 = obj.addEntityByDefaultID('VECTOR','',xn,1);
            v4 = obj.addEntityByDefaultID('VECTOR','',yn,1);

            l1 = obj.addEntityByDefaultID('LINE','',p1,v1);
            l2 = obj.addEntityByDefaultID('LINE','',p2,v2);
            l3 = obj.addEntityByDefaultID('LINE','',p3,v3);
            l4 = obj.addEntityByDefaultID('LINE','',p4,v4);

            v1 = obj.addEntityByDefaultID('VERTEX','',p1);
            v2 = obj.addEntityByDefaultID('VERTEX','',p2);
            v3 = obj.addEntityByDefaultID('VERTEX','',p3);
            v4 = obj.addEntityByDefaultID('VERTEX','',p4);

            e1 = obj.addEntityByDefaultID('EDGE_CURVE','',v1,v2,l1,1);
            e2 = obj.addEntityByDefaultID('EDGE_CURVE','',v2,v3,l2,1);
            e3 = obj.addEntityByDefaultID('EDGE_CURVE','',v3,v4,l3,1);
            e4 = obj.addEntityByDefaultID('EDGE_CURVE','',v4,v1,l4,1);

            e1 = obj.addEntityByDefaultID('ORIENTED_EDGE','',e1,1);
            e2 = obj.addEntityByDefaultID('ORIENTED_EDGE','',e2,1);
            e3 = obj.addEntityByDefaultID('ORIENTED_EDGE','',e3,1);
            e4 = obj.addEntityByDefaultID('ORIENTED_EDGE','',e4,1);

            loop1 = obj.addEntityByDefaultID('EDGE_LOOP','',[e1,e2,e3,e4]);
            b1 = obj.addEntityByDefaultID('FACE_OUTER_BOUND','',loop1,1);

            f1 = obj.addEntityByDefaultID('ADVANCED_FACE','',b1,xy_plane,1);

            % Closed shell containing the advanced face
            closedShell = obj.addEntityByDefaultID('CLOSED_SHELL', '', f1);

            % Manifold solid BREP referencing the closed shell
            manifoldSolid = obj.addEntityByDefaultID('MANIFOLD_SOLID_BREP', '', closedShell);

            % 3D geometric representation context
            geomContext = obj.addEntityByDefaultID('GEOMETRIC_REPRESENTATION_CONTEXT', '', 3);

            % Shape representation containing the manifold solid and context
            obj.addEntityByDefaultID('SHAPE_REPRESENTATION', '', manifoldSolid, geomContext);

        end

        function addCircle(obj, x0, y0, r)
            % Directions and axis placement (XY plane)
            origin = obj.addEntityByDefaultID('CARTESIAN_POINT','',[x0,y0,0]);
            x_dir = obj.addEntityByDefaultID('DIRECTION','',[1,0,0]);
            z_dir = obj.addEntityByDefaultID('DIRECTION','',[0,0,1]);
            axis_placement = obj.addEntityByDefaultID('AXIS2_PLACEMENT_3D','',origin, z_dir, x_dir);
            xy_plane = obj.addEntityByDefaultID('PLANE','',axis_placement);

            % Circle geometry
            circle = obj.addEntityByDefaultID('CIRCLE','',axis_placement, r);

            % Define 4 points on the circle (0°, 90°, 180°, 270°)
            pts = [
                x0 + r, y0, 0;
                x0, y0 + r, 0;
                x0 - r, y0, 0;
                x0, y0 - r, 0;
                ];
            verts = zeros(1,4);
            for i = 1:4
                cp = obj.addEntityByDefaultID('CARTESIAN_POINT','', pts(i,:) );
                verts(i) = obj.addEntityByDefaultID('VERTEX','',cp);
            end

            % Create 4 edges (quarter arcs) between the points using the circle
            edges = zeros(1,4);
            for i = 1:4
                start_v = verts(i);
                end_v = verts(mod(i,4) + 1); % next vertex, wrapping around
                edges(i) = obj.addEntityByDefaultID('EDGE_CURVE','',start_v, end_v, circle, true);
            end

            % Create oriented edges with positive orientation
            oriented_edges = zeros(1,4);
            for i = 1:4
                oriented_edges(i) = obj.addEntityByDefaultID('ORIENTED_EDGE','',edges(i), true);
            end

            % Create edge loop from oriented edges
            edge_loop = obj.addEntityByDefaultID('EDGE_LOOP','',oriented_edges);

            % Face outer bound with orientation true
            face_bound = obj.addEntityByDefaultID('FACE_OUTER_BOUND','',edge_loop, true);

            % Advanced face with the face boundary and the plane
            advanced_face = obj.addEntityByDefaultID('ADVANCED_FACE','',face_bound, xy_plane, true);

            % Closed shell containing this face
            closed_shell = obj.addEntityByDefaultID('CLOSED_SHELL', '', advanced_face);

            % Manifold solid BREP referencing the closed shell
            manifold_solid = obj.addEntityByDefaultID('MANIFOLD_SOLID_BREP', '', closed_shell);

            % 3D geometric representation context
            geom_context = obj.addEntityByDefaultID('GEOMETRIC_REPRESENTATION_CONTEXT', '', 3);

            % Shape representation containing the solid and context
            obj.addEntityByDefaultID('SHAPE_REPRESENTATION', '', manifold_solid, geom_context);
        end

        function addPolygon(obj, xCoords, yCoords)
            % Validate input lengths
            n = length(xCoords);
            if n < 3 || n ~= length(yCoords)
                error('xCoords and yCoords must have same length >= 3');
            end

            % Define directions and plane (XY plane)
            origin = obj.addEntityByDefaultID('CARTESIAN_POINT','',[0,0,0]);
            x_dir = obj.addEntityByDefaultID('DIRECTION','',[1,0,0]);
            y_dir = obj.addEntityByDefaultID('DIRECTION','',[0,1,0]);
            z_dir = obj.addEntityByDefaultID('DIRECTION','',[0,0,1]);
            axis_placement = obj.addEntityByDefaultID('AXIS2_PLACEMENT_3D','',origin, z_dir, x_dir);
            xy_plane = obj.addEntityByDefaultID('PLANE','',axis_placement);

            % Create Cartesian points for polygon vertices
            pts = zeros(1,n);
            verts = zeros(1,n);
            for i = 1:n
                pts(i) = obj.addEntityByDefaultID('CARTESIAN_POINT','',[xCoords(i), yCoords(i), 0]);
                verts(i) = obj.addEntityByDefaultID('VERTEX','',pts(i));
            end

            % Create vectors for each edge (from vertex i to i+1)
            % For edges, direction vector = next point - current point
            % We'll create a VECTOR and a LINE for each edge
            edges = zeros(1,n);
            oriented_edges = zeros(1,n);

            for i = 1:n
                idxNext = mod(i, n) + 1; % next vertex index, wrapping around

                % Direction vector = next point - current point
                dx = xCoords(idxNext) - xCoords(i);
                dy = yCoords(idxNext) - yCoords(i);

                dir_vec = obj.addEntityByDefaultID('DIRECTION','',[dx, dy, 0]);
                vec = obj.addEntityByDefaultID('VECTOR','',dir_vec, sqrt(dx^2 + dy^2));

                % Create LINE entity for this edge
                edges(i) = obj.addEntityByDefaultID('LINE','',pts(i), vec);

                % Create EDGE_CURVE entity
                edge_curve = obj.addEntityByDefaultID('EDGE_CURVE','',verts(i), verts(idxNext), edges(i), true);

                % Create ORIENTED_EDGE entity with orientation true
                oriented_edges(i) = obj.addEntityByDefaultID('ORIENTED_EDGE','',edge_curve, true);
            end

            % Create EDGE_LOOP from all oriented edges
            edge_loop = obj.addEntityByDefaultID('EDGE_LOOP','',oriented_edges);

            % FACE_OUTER_BOUND with orientation true
            face_bound = obj.addEntityByDefaultID('FACE_OUTER_BOUND','',edge_loop, true);

            % ADVANCED_FACE from the boundary and plane
            advanced_face = obj.addEntityByDefaultID('ADVANCED_FACE','',face_bound, xy_plane, true);

            % Closed shell containing the face
            closed_shell = obj.addEntityByDefaultID('CLOSED_SHELL', '', advanced_face);

            % Manifold solid BREP referencing the closed shell
            manifold_solid = obj.addEntityByDefaultID('MANIFOLD_SOLID_BREP', '', closed_shell);

            % 3D geometric representation context
            geom_context = obj.addEntityByDefaultID('GEOMETRIC_REPRESENTATION_CONTEXT', '', 3);

            % Shape representation containing the solid and context
            obj.addEntityByDefaultID('SHAPE_REPRESENTATION', '', manifold_solid, geom_context);
        end
        
        function addBox(obj, x0, y0, z0, w, h, d)
            
            % Directions used as references
            xp = obj.addEntityByDefaultID('DIRECTION','',[1,0,0]);
            yp = obj.addEntityByDefaultID('DIRECTION','',[0,1,0]);
            zp = obj.addEntityByDefaultID('DIRECTION','',[0,0,1]);
            xn = obj.addEntityByDefaultID('DIRECTION','',[-1,0,0]);
            yn = obj.addEntityByDefaultID('DIRECTION','',[0,-1,0]);
            zn = obj.addEntityByDefaultID('DIRECTION','',[0,0,-1]);

            % Helper function to create plane entity with placement
            function plane_id = createPlane(origin_xyz, normal_vec, ref_vec)
                origin = obj.addEntityByDefaultID('CARTESIAN_POINT','',origin_xyz);
                normal = obj.addEntityByDefaultID('DIRECTION','',normal_vec);
                ref = obj.addEntityByDefaultID('DIRECTION','',ref_vec);
                axis_placement = obj.addEntityByDefaultID('AXIS2_PLACEMENT_3D','',origin, normal, ref);
                plane_id = obj.addEntityByDefaultID('PLANE','',axis_placement);
            end

            % Define the 8 box corner points
            pts = zeros(1,8);
            pts(1) = obj.addEntityByDefaultID('CARTESIAN_POINT','',[x0      , y0      , z0     ]);
            pts(2) = obj.addEntityByDefaultID('CARTESIAN_POINT','',[x0 + w  , y0      , z0     ]);
            pts(3) = obj.addEntityByDefaultID('CARTESIAN_POINT','',[x0 + w  , y0 + h  , z0     ]);
            pts(4) = obj.addEntityByDefaultID('CARTESIAN_POINT','',[x0      , y0 + h  , z0     ]);
            pts(5) = obj.addEntityByDefaultID('CARTESIAN_POINT','',[x0      , y0      , z0 + d ]);
            pts(6) = obj.addEntityByDefaultID('CARTESIAN_POINT','',[x0 + w  , y0      , z0 + d ]);
            pts(7) = obj.addEntityByDefaultID('CARTESIAN_POINT','',[x0 + w  , y0 + h  , z0 + d ]);
            pts(8) = obj.addEntityByDefaultID('CARTESIAN_POINT','',[x0      , y0 + h  , z0 + d ]);

            % Create vectors for edges (from directions)
            vectors = zeros(1,12);
            % Each edge vector is based on direction and length between points
            % You can assign them similarly to rectangle code by direction and length 1

            % We'll create 12 edge vectors (along box edges)
            vectors(1) = obj.addEntityByDefaultID('VECTOR','',xp, w);
            vectors(2) = obj.addEntityByDefaultID('VECTOR','',yp, h);
            vectors(3) = obj.addEntityByDefaultID('VECTOR','',zp, d);
            vectors(4) = obj.addEntityByDefaultID('VECTOR','',xn, w);
            vectors(5) = obj.addEntityByDefaultID('VECTOR','',yn, h);
            vectors(6) = obj.addEntityByDefaultID('VECTOR','',zn, d);

            % Create 12 lines (edges) connecting points and vectors
            lines = zeros(1,12);
            lines(1)  = obj.addEntityByDefaultID('LINE','',pts(1),vectors(1));
            lines(2)  = obj.addEntityByDefaultID('LINE','',pts(2),vectors(2));
            lines(3)  = obj.addEntityByDefaultID('LINE','',pts(3),vectors(4));
            lines(4)  = obj.addEntityByDefaultID('LINE','',pts(4),vectors(5));
            lines(5)  = obj.addEntityByDefaultID('LINE','',pts(5),vectors(1));
            lines(6)  = obj.addEntityByDefaultID('LINE','',pts(6),vectors(2));
            lines(7)  = obj.addEntityByDefaultID('LINE','',pts(7),vectors(4));
            lines(8)  = obj.addEntityByDefaultID('LINE','',pts(8),vectors(5));
            lines(9)  = obj.addEntityByDefaultID('LINE','',pts(1),vectors(3));
            lines(10) = obj.addEntityByDefaultID('LINE','',pts(2),vectors(3));
            lines(11) = obj.addEntityByDefaultID('LINE','',pts(3),vectors(3));
            lines(12) = obj.addEntityByDefaultID('LINE','',pts(4),vectors(3));

            % Create vertices for each point
            vertices = zeros(1,8);
            for i=1:8
                vertices(i) = obj.addEntityByDefaultID('VERTEX','',pts(i));
            end

            % Create 12 EDGE_CURVE entities linking vertices and lines
            edge_curves = zeros(1,12);
            edge_curves(1)  = obj.addEntityByDefaultID('EDGE_CURVE','',vertices(1), vertices(2), lines(1), true);
            edge_curves(2)  = obj.addEntityByDefaultID('EDGE_CURVE','',vertices(2), vertices(3), lines(2), true);
            edge_curves(3)  = obj.addEntityByDefaultID('EDGE_CURVE','',vertices(3), vertices(4), lines(3), true);
            edge_curves(4)  = obj.addEntityByDefaultID('EDGE_CURVE','',vertices(4), vertices(1), lines(4), true);
            edge_curves(5)  = obj.addEntityByDefaultID('EDGE_CURVE','',vertices(5), vertices(6), lines(5), true);
            edge_curves(6)  = obj.addEntityByDefaultID('EDGE_CURVE','',vertices(6), vertices(7), lines(6), true);
            edge_curves(7)  = obj.addEntityByDefaultID('EDGE_CURVE','',vertices(7), vertices(8), lines(7), true);
            edge_curves(8)  = obj.addEntityByDefaultID('EDGE_CURVE','',vertices(8), vertices(5), lines(8), true);
            edge_curves(9)  = obj.addEntityByDefaultID('EDGE_CURVE','',vertices(1), vertices(5), lines(9), true);
            edge_curves(10) = obj.addEntityByDefaultID('EDGE_CURVE','',vertices(2), vertices(6), lines(10), true);
            edge_curves(11) = obj.addEntityByDefaultID('EDGE_CURVE','',vertices(3), vertices(7), lines(11), true);
            edge_curves(12) = obj.addEntityByDefaultID('EDGE_CURVE','',vertices(4), vertices(8), lines(12), true);

            % Create ORIENTED_EDGE entities for all edges (all with true orientation)
            oriented_edges = zeros(1,12);
            for i=1:12
                oriented_edges(i) = obj.addEntityByDefaultID('ORIENTED_EDGE','',edge_curves(i), true);
            end

            % Helper to create edge loop from oriented edges
            function loop_id = createLoop(edge_ids)
                loop_id = obj.addEntityByDefaultID('EDGE_LOOP','',edge_ids);
            end

            % Create planes for each face:

            % Bottom face plane (XY, z=z0)
            planes(1) = createPlane([x0, y0, z0], [0,0,-1], [1,0,0]);
            % Top face plane (XY, z=z0+d)
            planes(2) = createPlane([x0, y0, z0 + d], [0,0,1], [1,0,0]);
            % Front face plane (XZ, y=y0)
            planes(3) = createPlane([x0, y0, z0], [0,-1,0], [1,0,0]);
            % Back face plane (XZ, y=y0+h)
            planes(4) = createPlane([x0, y0 + h, z0], [0,1,0], [1,0,0]);
            % Left face plane (YZ, x=x0)
            planes(5) = createPlane([x0, y0, z0], [-1,0,0], [0,0,1]);
            % Right face plane (YZ, x=x0+w)
            planes(6) = createPlane([x0 + w, y0, z0], [1,0,0], [0,0,1]);

            % Create edge loops for each face (by grouping corresponding oriented edges)
            % Indices for edges per face:
            % Bottom: edges 1,2,3,4 (loop bottom)
            loop_bottom = createLoop(oriented_edges([1,2,3,4]));
            % Top: edges 5,6,7,8
            loop_top = createLoop(oriented_edges([5,6,7,8]));
            % Front: edges 1,10,5,9 (edges along front face)
            loop_front = createLoop(oriented_edges([1,10,5,9]));
            % Back: edges 3,12,7,11
            loop_back = createLoop(oriented_edges([3,12,7,11]));
            % Left: edges 4,9,8,12
            loop_left = createLoop(oriented_edges([4,9,8,12]));
            % Right: edges 2,11,6,10
            loop_right = createLoop(oriented_edges([2,11,6,10]));

            % Create face bounds from loops
            face_bounds = zeros(1,6);
            face_bounds(1) = obj.addEntityByDefaultID('FACE_OUTER_BOUND','',loop_bottom, true);
            face_bounds(2) = obj.addEntityByDefaultID('FACE_OUTER_BOUND','',loop_top, true);
            face_bounds(3) = obj.addEntityByDefaultID('FACE_OUTER_BOUND','',loop_front, true);
            face_bounds(4) = obj.addEntityByDefaultID('FACE_OUTER_BOUND','',loop_back, true);
            face_bounds(5) = obj.addEntityByDefaultID('FACE_OUTER_BOUND','',loop_left, true);
            face_bounds(6) = obj.addEntityByDefaultID('FACE_OUTER_BOUND','',loop_right, true);

            % Create advanced faces, assigning each face its corresponding plane
            advanced_faces = zeros(1,6);
            for i=1:6
                advanced_faces(i) = obj.addEntityByDefaultID('ADVANCED_FACE','',face_bounds(i), planes(i), true);
            end

            % Create CLOSED_SHELL from advanced faces
            closed_shell = obj.addEntityByDefaultID('CLOSED_SHELL','',advanced_faces);

            % Create MANIFOLD_SOLID_BREP referencing closed shell
            manifold_solid = obj.addEntityByDefaultID('MANIFOLD_SOLID_BREP','',closed_shell);

            % Create GEOMETRIC_REPRESENTATION_CONTEXT for 3D
            geom_context = obj.addEntityByDefaultID('GEOMETRIC_REPRESENTATION_CONTEXT','',3);

            % Finally, create SHAPE_REPRESENTATION containing the manifold solid and geom context
            obj.addEntityByDefaultID('SHAPE_REPRESENTATION','',manifold_solid, geom_context);

        end

        function writeSTEP(obj, fileName)
            fid = fopen(fileName, 'w');
            if fid == -1
                error('Could not open file for writing: %s', fileName);
            end

            % Write the STEP header
            fprintf(fid, 'ISO-10303-21;\n');
            fprintf(fid, 'HEADER;\n');
            % %             fprintf(fid, "FILE_DESCRIPTION(('MATLAB Generated STEP'), '1');\n");
            % %             fprintf(fid, "FILE_NAME('%s', '%s');\n", fileName, datestr(now, 'yyyy-mm-ddTHH:MM:SS'));
            %             fprintf(fid, "FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\n");
            fprintf(fid, 'ENDSEC;\n');

            % Write the DATA section
            fprintf(fid, 'DATA;\n');
            for i = 1:numel(obj.entities)
                entity = obj.entities{i};
                if ismethod(entity, 'toSTEP')
                    line = entity.toSTEP();  % Convert to STEP text
                    fprintf(fid, '%s\n', line);
                else
                    warning('Entity ID %d does not implement toSTEP method.', entity.id);
                end
            end
            fprintf(fid, 'ENDSEC;\n');
            fprintf(fid, 'END-ISO-10303-21;\n');

            fclose(fid);
            fprintf('STEP file written to %s\n', fileName);
        end

    end
end