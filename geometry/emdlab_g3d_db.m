% developer: https://ComProgExpert.com
% data base class for 3D geometries

classdef emdlab_g3d_db < handle

    properties

        % points
        points (1,:) emdlab_g3d_point;

        % edges
        edges (1,:) emdlab_g3d_edge;

        % loops
        loops (1,:) emdlab_g3d_loop;

        % faces
        faces (1,:) emdlab_g3d_face;

    end

    methods
        %% constructor and destructor
        function obj = emdlab_g3d_db()
        end

        %% point methods
        % adding a new point to data base
        % this function returns point index and point handle
        function varargout = addPoint(obj, varargin)

            % check the varargin type
            if numel(varargin) == 2
                x = varargin{1};
                y = varargin{2};
            elseif numel(varargin) == 1
                if isa(varargin{1},'emdlab_g2d_point')
                    x = varargin{1}.x;
                    y = varargin{1}.y;
                else
                    throw(MException('', 'Wrong input type, it must be <emdlab_g2d_point> type.'));
                end
            else
                throw(MException('', 'Wrong number of input arguments.'));
            end

            % check for existance of already defined point in the same location
            for i = 1:numel(obj.points)

                if norm(obj.points(i).getVector() - [x,y]) < 1e-6

                    pointHandle = obj.points(i);
                    if nargout == 1
                        varargout{1} = i;
                    elseif nargout == 2
                        varargout{1} = i;
                        varargout{2} = pointHandle;
                    elseif nargout > 2
                        error('The number of output arguments is too high.');
                    end
                    return;

                end

            end

            % get an instance of point class
            pointHandle = emdlab_g2d_point(x,y);
            obj.points(end+1) = pointHandle;

            % generate point tag
            pointIndex = numel(obj.points);
            pointHandle.tag = ['p', num2str(pointIndex)];

            if nargout == 1
                varargout{1} = pointIndex;
            elseif nargout == 2
                varargout{1} = pointIndex;
                varargout{2} = pointHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        %% edge methods
        % adding a new segment to data base
        % this function returns edge index and edge handle
        function varargout = addSegment(obj, p0Index, p1Index)

            % get an instance of edge class
            edgeHandle = emdlab_g2d_edge;
            % set pointer class to segment
            edgeHandle.ptr = emdlab_g2d_segment(obj.points(p0Index),obj.points(p1Index));
            obj.edges(end+1) = edgeHandle;
                       
            % generate edge tag
            edgeIndex = numel(obj.edges);
            edgeHandle.tag = ['e', num2str(edgeIndex)];

            if nargout == 1
                varargout{1} = edgeIndex;
            elseif nargout == 2
                varargout{1} = edgeIndex;
                varargout{2} = edgeHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % add a new segment by direct coordinates passing
        function varargout = addSegmentByCoordinates(obj, x1, y1, x2, y2)

            p1Index = obj.addPoint(x1,y1);
            p2Index = obj.addPoint(x2,y2);

            if nargout == 0
                obj.addSegment(p1Index,p2Index);
            elseif nargout == 1
                varargout{1} = obj.addSegment(p1Index,p2Index);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSegment(p1Index,p2Index);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % adding a new arc to data base
        function varargout = addArc(obj, p0Index, p1Index, p2Index, direction)

            % get an instance of edge class
            edgeHandle = emdlab_g2d_edge;
            % set pointer class to arc
            edgeHandle.ptr = emdlab_g2d_arc(obj.points(p0Index),obj.points(p1Index),obj.points(p2Index), direction);
            obj.edges(end+1) = edgeHandle;

            % generate edge tag
            edgeIndex = numel(obj.edges);
            edgeHandle.tag = ['e', num2str(edgeIndex)];

            if nargout == 1
                varargout{1} = edgeIndex;
            elseif nargout == 2
                varargout{1} = edgeIndex;
                varargout{2} = edgeHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % get edge handle
        function edgeHandle = getEdgeHandleByIndex(obj, eIndex)

            % check for existance of already defined edge in data base
            if eIndex <= numel(obj.edges)
                edgeHandle = obj.edges(eIndex);
            else
                error('Edge was not found.');
            end

        end

        function edgeHandle = getEdgeHandleByTag(obj, eTag)

            % check for existance of already defined edge in data base
            for i = 1:numel(obj.edges)

                if strcmp(obj.edges(i).tag,eTag)

                    edgeHandle = obj.edges(i);
                    return;

                end

            end

            error('Edge was not found.');

        end

        %% loop methods
        % adding a new loop to data base
        % this function returns loop index and loop handle
        function varargout = addLoop(obj, varargin)

            % get a loop class instance
            loopHandle = emdlab_g2d_loop;
            obj.loops(end+1) = loopHandle;

            % loop to assign edges and directions
            for i = 1:2:numel(varargin)

                loopHandle.addEdge(obj.edges(varargin{i}).ptr, varargin{i+1});

            end

            % generate loop tag
            loopIndex = numel(obj.loops);
            loopHandle.tag = ['l', num2str(loopIndex)];

            if nargout == 1
                varargout{1} = loopIndex;
            elseif nargout == 2
                varargout{1} = loopIndex;
                varargout{2} = loopHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        %% face methods
        % adding a new face to data base
        % this function returns the face handle
        function varargout = addFace(obj, faceName, varargin)

            % get face class instance
            faceHandle = emdlab_g2d_face;
            obj.faces(end+1) = faceHandle;

            for i = 1:numel(varargin)

                obj.faces(end).addLoop(obj.loops(varargin{i}));

            end
            faceHandle.tag = faceName;

            if nargout == 1
                varargout = faceHandle;
            elseif nargout > 1
                error('The number of output arguments is too hight');
            end

        end

        %% mesh generation methods
        % generate triangular mesh for geometry
        function m = generateMesh(obj, meshGenerator)

            % default mesh generator
            if nargin<2
                meshGenerator = 'mm';
            end

            % get an instance of mesh data base
            m = emdlab_m2d_tmdb;

            % add mesh zones
            for i = 1:numel(obj.faces)
                m.addMeshZone(obj.faces(i).tag, obj.faces(i).getMesh(meshGenerator));
            end

        end

        %% visualization methos
        % show the geometry sketch
        function varargout = showSketch(obj, showTags)

            if nargin<2
                showTags = true;
            end

            f = figure('NumberTitle', 'on', 'WindowState', 'maximized', 'name', 'EMDLAB Geometry Visualization', 'color', [0.9,0.9,0.9]);
            hold all;

            % plot points: Np = the number of points
            Np = numel(obj.points);
            p = zeros(Np,2);
            for i = 1:Np
                p(i,1) = obj.points(i).x;
                p(i,2) = obj.points(i).y;
            end
            if ~isempty(p)
                pointTags = cell(1,Np);
                for i = 1:Np
                    pointTags{i} = obj.points(i).tag;
                end
                if showTags
                    text(p(:,1), p(:,2), pointTags, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'BackgroundColor', 'y');
                end
                plot(p(:,1), p(:,2), 's', 'LineWidth', 1.5);
                p_min = min(p,[],1);
                p_max = max(p,[],1);
                plot([min(p_min(1),0),max(p_max(1),0)], [0,0], '--', 'color', [0.2,0.2,0.2]);
                plot([0,0], [min(p_min(2),0),max(p_max(2),0)], '--', 'color', [0.2,0.2,0.2]);
            end

            % plot edges: Ne = the number of edges
            Ne = numel(obj.edges);
            v = cell(Ne,1);
            cl = cell(Ne,1);
            c = zeros(Ne,2);
            for i = 1:Ne
                v{i} = obj.edges(i).ptr.getMeshNodes;
                cl{i} = (1:size(v{i},1)-1)';
                cl{i} = [cl{i},cl{i}+1];
                c(i,:) = obj.edges(i).ptr.getCenter;
            end
            Index = 0;
            for i = 2:Ne
                Index = Index + size(v{i-1},1);
                cl{i} = cl{i} + Index;
            end

            v = cell2mat(v);
            cl = cell2mat(cl);

            if ~isempty(v)
                patch('faces', cl, 'vertices', v, 'edgecolor', 'b', 'linewidth',1.2);
                plot(v(:,1), v(:,2), '.', 'color', 'k');
                if showTags
                    edgeTags = cell(1,Ne);
                    for i = 1:Ne
                        edgeTags{i} = obj.edges(i).tag;
                    end
                    text(c(:,1), c(:,2), edgeTags, 'BackgroundColor', 'w', ...
                        'HorizontalAlignment','center','VerticalAlignment','middle');
                end
            end

            set(gca, 'clipping', 'off');

            axis off equal
            zoom on
            grid on

            if nargout == 1
                varargout{1} = f;
            end

        end

        % adding primitive loops
        % this function returns loop index and loop handle
        function varargout = addRectangleLoop(obj, x0, y0, w, h)

            p1Index = obj.addPoint(x0,y0);
            p2Index = obj.addPoint(x0+w,y0);
            p3Index = obj.addPoint(x0+w,y0+h);
            p4Index = obj.addPoint(x0,y0+h);

            e1Index = obj.addSegment(p1Index, p2Index);
            e2Index = obj.addSegment(p2Index, p3Index);
            e3Index = obj.addSegment(p3Index, p4Index);
            e4Index = obj.addSegment(p4Index, p1Index);

            if nargout == 0
                obj.addLoop(e1Index, 1, e2Index, 1, e3Index, 1, e4Index, 1);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, 1, e2Index, 1, e3Index, 1, e4Index, 1);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, 1, e2Index, 1, e3Index, 1, e4Index, 1);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % this function returns loop index and loop handle
        function varargout = addCircleLoop(obj, x0, y0, r)

            p1Index = obj.addPoint(x0,y0);
            p2Index = obj.addPoint(x0+r,y0);
            p3Index = obj.addPoint(x0-r,y0);

            e1Index = obj.addArc(p1Index, p2Index, p3Index, 1);
            e2Indexe = obj.addArc(p1Index, p3Index, p2Index, 1);

            if nargout == 0
                obj.addLoop(e1Index, 1, e2Indexe, 1);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, 1, e2Indexe, 1);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, 1, e2Indexe, 1);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % this function returns loop index and loop handle
        function varargout = addClosedPolylineLoop(obj, x, y)

            Nx = length(x);
            Ny = length(y);

            if Nx~=Ny
                error('Number of x and y coordinates must be the same.');
            end

            p_indices = zeros(1,Nx);
            for i = 1:length(x)
                p_indices(i) = obj.addPoint(x(i),y(i));
            end

            p_indices(end+1) = p_indices(1);
            e_indices = cell(1,2*Nx);
            for i = 1:Nx
                e_indices{2*i-1} = obj.addSegment(p_indices(i),p_indices(i+1));
                e_indices{2*i} = 1;
            end

            if nargout == 0
                obj.addLoop(e_indices{:});
            elseif nargout == 1
                varargout{1} = obj.addLoop(e_indices{:});
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e_indices{:});
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % mesh max length
        function setMeshMaxLength(obj, mLength)

            for i = 1:numel(obj.edges)
                obj.edges(i).ptr.setMaxLength(mLength);
            end

        end

        function setMeshLengthByRadialFunction(obj, fHandle)

            for i = 1:numel(obj.edges)

                if isa(obj.edges(i).ptr, 'emdlab_g2d_segment')
                    obj.edges(i).ptr.setL1L2(fHandle(obj.edges(i).ptr.p0.norm()), fHandle(obj.edges(i).ptr.p1.norm()))
                else
                    obj.edges(i).ptr.setMaxLength(fHandle(obj.edges(i).ptr.p1.norm()));
                end

            end

        end

    end

end