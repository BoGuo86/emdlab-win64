% EMDLAB: Electrical Machines Design Laboratory
% data base class for 2d geometries

classdef emdlab_g2d_db < handle

    properties

        % points
        points (1,:) emdlab_g2d_point;

        % edges
        edges (1,:) emdlab_g2d_edge;

        % loops
        loops (1,:) emdlab_g2d_loop;

        % faces
        faces (1,:) emdlab_g2d_face;

        % python path
        pyPath = "";

    end

    methods
        %% constructor and destructor
        function obj = emdlab_g2d_db()

            % set python path
            fid = fopen('c:geometry\pyPath.txt', 'r');
            if fid<0
                warning('Cannot set Python path. The file <pyPath.txt> does not exist.');
            else
                obj.pyPath = string(fgetl(fid));
                obj.pyPath = replace(obj.pyPath, '\', '\\');
                fclose(fid);
            end

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

        function removePoint(obj, pIndex)

            % remove connected edges to this point firt
            for i = 1:numel(pIndex.tags)
            end
            obj.points(pIndex) = [];
        end

        function pointHandle = getPointHandleByTag(obj, pTag)

            % check for existance of already defined point in data base
            for i = 1:numel(obj.points)

                if strcmp(obj.points(i).tag,pTag)

                    pointHandle = obj.points(i);
                    return;

                end

            end

            error('Point was not found.');

        end

        function pointIndex = getPointIndexByTag(obj, pTag)

            % check for existance of already defined point in data base
            for i = 1:numel(obj.points)

                if strcmp(obj.points(i).tag,pTag)

                    pointIndex = i;
                    return;

                end

            end

            error('Point was not found.');

        end

        function pointIndex = getPointIndexByCoordinates(obj, x, y)

            if ~numel(obj.points)
                error('There is no defined point.');
            end

            % this function returns index of closed point to x and y coordinates
            minDistance = inf;
            for i = 1:numel(obj.points)

                distance = norm(obj.points(i).getVector - [x,y]);
                if distance < minDistance
                    pointIndex = i;
                    minDistance = distance;
                end

            end

        end

        function alignPointsAlongYAxis(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).x = obj.points(varargin{1}).x;
            end

        end

        function alignPointsAlongXAxis(obj, varargin)

            for i = 2:numel(varargin)
                obj.points(varargin{i}).y = obj.points(varargin{1}).y;
            end

        end

        function alignPointsAlongRAxis(obj, varargin)

            u_ref = obj.points(varargin{1}).getUnitVector;
            for i = 2:numel(varargin)
                r_i = obj.points(varargin{i}).getDistanceFromOrigin;
                obj.points(varargin{i}).setCoordinates(r_i*u_ref(1), r_i*u_ref(2));
            end

        end

        function alignPointsAlongTAxis(obj, varargin)

            r_ref = obj.points(varargin{1}).getDistanceFromOrigin;
            for i = 2:numel(varargin)
                u_i = obj.points(varargin{i}).getUnitVector;
                obj.points(varargin{i}).setCoordinates(r_ref*u_i(1), r_ref*u_i(2));
            end

        end

        %% edge methods
        % adding a new segment to data base
        % this function returns edge index and edge handle
        function varargout = addSegment(obj, p0Index, p1Index)

            % check for char inputs
            if ischar(p0Index)
                p0Index = obj.getPointIndexByTag(p0Index);
            end

            if ischar(p1Index)
                p1Index = obj.getPointIndexByTag(p1Index);
            end

            % get an instance of edge class
            edgeHandle = emdlab_g2d_edge;
            % set pointer class to segment
            edgeHandle.ptr = emdlab_g2d_segment(obj.points(p0Index),obj.points(p1Index));
            obj.edges(end+1) = edgeHandle;

            % generate edge tag
            edgeIndex = numel(obj.edges);
            edgeHandle.tag = ['e', num2str(edgeIndex)];
            edgeHandle.ptr.tag = edgeHandle.tag;

            % add edge tag to connected points
            obj.points(p0Index).tags(end+1) = edgeHandle.tag;
            obj.points(p1Index).tags(end+1) = edgeHandle.tag;

            if nargout == 1
                varargout{1} = edgeIndex;
            elseif nargout == 2
                varargout{1} = edgeIndex;
                varargout{2} = edgeHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % adding a new spline to data base
        % this function returns edge index and edge handle
        function varargout = addSpline(obj, ptsIndex)

            % number of points
            Np = numel(ptsIndex);
            pts = repmat(emdlab_g2d_point,1,Np);

            if isa(ptsIndex, 'cell')

                for i = 1:Np
                    if ischar(ptsIndex{i})
                        pts(i) = obj.points(obj.getPointIndexByTag(ptsIndex{i}));
                    else
                        pts(i) = obj.points(ptsIndex{i});
                    end
                end

            else

                for i = 1:Np
                    pts(i) = obj.points(ptsIndex(i));
                end

            end

            % get an instance of edge class
            edgeHandle = emdlab_g2d_edge;

            % set pointer class to segment
            edgeHandle.ptr = emdlab_g2d_spline(pts);
            obj.edges(end+1) = edgeHandle;

            % generate edge tag
            edgeIndex = numel(obj.edges);
            edgeHandle.tag = ['e', num2str(edgeIndex)];
            edgeHandle.ptr.tag = edgeHandle.tag;

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

        % add a new segment by direct coordinates passing
        function varargout = addSplineByCoordinates(obj, x, y)

            % check feasibility
            if length(x) ~= length(y)
                error('x and y vectors must have the same length.');
            end

            Nx = length(x);
            ptsIndex = zeros(1,Nx);
            for i = 1:Nx
                ptsIndex(i) = obj.addPoint(x(i),y(i));
            end

            if nargout == 0
                obj.addSpline(ptsIndex);
            elseif nargout == 1
                varargout{1} = obj.addSpline(ptsIndex);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSpline(ptsIndex);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % adding a new arc to data base
        function varargout = addArc(obj, p0Index, p1Index, p2Index, direction)

            % check for char inputs
            if ischar(p0Index)
                p0Index = obj.getPointIndexByTag(p0Index);
            end

            if ischar(p1Index)
                p1Index = obj.getPointIndexByTag(p1Index);
            end

            if ischar(p2Index)
                p2Index = obj.getPointIndexByTag(p2Index);
            end

            % get an instance of edge class
            edgeHandle = emdlab_g2d_edge;
            % set pointer class to arc
            edgeHandle.ptr = emdlab_g2d_arc(obj.points(p0Index),obj.points(p1Index),obj.points(p2Index), direction);
            obj.edges(end+1) = edgeHandle;

            % generate edge tag
            edgeIndex = numel(obj.edges);
            edgeHandle.tag = ['e', num2str(edgeIndex)];
            edgeHandle.ptr.tag = edgeHandle.tag;

            % add edge to point tags
            obj.points(p0Index).tags(end+1) = edgeHandle.tag;
            obj.points(p1Index).tags(end+1) = edgeHandle.tag;
            obj.points(p2Index).tags(end+1) = edgeHandle.tag;

            if nargout == 1
                varargout{1} = edgeIndex;
            elseif nargout == 2
                varargout{1} = edgeIndex;
                varargout{2} = edgeHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addArcByCoordinates(obj, x1, y1, x2, y2, x3, y3, direction)

            p1Index = obj.addPoint(x1,y1);
            p2Index = obj.addPoint(x2,y2);
            p3Index = obj.addPoint(x3,y3);

            if nargout == 0
                obj.addArc(p1Index, p2Index, p3Index, direction);
            elseif nargout == 1
                varargout{1} = obj.addArc(p1Index, p2Index, p3Index, direction);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addArc(p1Index, p2Index, p3Index, direction);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addArcByCoordinatesCPA(obj, x1, y1, x2, y2, arcAngle)

            p1Index = obj.addPoint(x1,y1);
            p2Index = obj.addPoint(x2,y2);
            [x3,y3] = emdlab_g2d_rotatePointsXY(x2,y2,arcAngle,x1,y1);
            p3Index = obj.addPoint(x3,y3);
            direction = arcAngle > 0;

            if nargout == 0
                obj.addArc(p1Index, p2Index, p3Index, direction);
            elseif nargout == 1
                varargout{1} = obj.addArc(p1Index, p2Index, p3Index, direction);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addArc(p1Index, p2Index, p3Index, direction);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        % get edge handle
        function edgeHandle = getEdgeHandleByIndex(obj, eIndex)

            % check for existance of already defined edge in data base
            if eIndex <= numel(obj.edges)
                edgeHandle = obj.edges(eIndex).ptr;
            else
                error('Edge was not found.');
            end

        end

        function edgeHandle = getEdgeHandleByTag(obj, eTag)

            % check for existance of already defined edge in data base
            for i = 1:numel(obj.edges)

                if strcmp(obj.edges(i).tag,eTag)

                    edgeHandle = obj.edges(i).ptr;
                    return;

                end

            end

            error('Edge was not found.');

        end

        function edgeIndex = getEdgeIndexByTag(obj, eTag)

            % check for existance of already defined edge in data base
            for i = 1:numel(obj.edges)

                if strcmp(obj.edges(i).tag,eTag)

                    edgeIndex = i;
                    return;

                end

            end

            error('Edge was not found.');

        end

        function splitSegment(obj, eIndex)

            edgeHandle = obj.edges(eIndex).ptr;
            tmp = edgeHandle.getCenter;
            p = obj.addPoint(tmp(1),tmp(2));
            p2 = edgeHandle.p1;
            edgeHandle.p1 = obj.points(p);
            obj.addSegment(p,obj.addPoint(p2));

        end

        function varargout = extendSegmentBySegment(obj, eIndex, extAngle, extAmplitude, seIndex)

            % set default start/end index as end
            if nargin < 5, seIndex = 1; end

            % get segment edge handle
            edgeHandle = obj.edges(eIndex).ptr;
            u = edgeHandle.getUnitVector;
            u = emdlab_g2d_rotatePoints(u,extAngle);

            % check start or index: 0 when we exten inward, 0: outward extension
            if seIndex == 0
                v2 = edgeHandle.p0.getVector;
                v1 = v2 + extAmplitude * u;
            else
                v1 = edgeHandle.p1.getVector;
                v2 = v1 + extAmplitude * u;
            end

            if nargout == 0
                obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 1
                varargout{1} = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = extendSegmentByTangentArc(obj, eIndex, extRadius, extAngle, seIndex)

            % set default start/end index as end
            if nargin < 5, seIndex = 1; end

            edgeHandle = obj.edges(eIndex).ptr;
            u = edgeHandle.getUnitVector;
            u = emdlab_g2d_rotatePoints(u,pi/2);

            % check start or index: 0 when we exten inward, 0: outward extension
            if seIndex == 0
                v2 = edgeHandle.p0.getVector;
                c = v2 + extRadius * u;
                v1 = emdlab_g2d_rotatePoints(v2, -sign(extRadius)*extAngle, c(1), c(2));
            else
                v1 = edgeHandle.p1.getVector;
                c = v1 + extRadius * u;
                v2 = emdlab_g2d_rotatePoints(v1, sign(extRadius)*extAngle, c(1), c(2));
            end

            if nargout == 0
                obj.addArcByCoordinates(c(1),c(2),v1(1),v1(2),v2(1),v2(2),extRadius>0);
            elseif nargout == 1
                varargout{1} = obj.addArcByCoordinates(c(1),c(2),v1(1),v1(2),v2(1),v2(2),extRadius>0);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addArcByCoordinates(c(1),c(2),v1(1),v1(2),v2(1),v2(2),extRadius>0);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = extendSegmentByArc(obj, eIndex, xc, yc, extAngle, seIndex)

            % set default start/end index as end
            if nargin < 6, seIndex = 1; end

            % get edge handle
            edgeHandle = obj.edges(eIndex).ptr;

            % check start or end point of the segment for extension
            if seIndex == 0
                v1 = edgeHandle.p0.getVector;
            else
                v1 = edgeHandle.p1.getVector;
            end

            v2 = emdlab_g2d_rotatePoints(v1, extAngle, xc, yc);
            if nargout == 0
                obj.addArcByCoordinates(xc,yc,v1(1),v1(2),v2(1),v2(2),extAngle>0);
            elseif nargout == 1
                varargout{1} = obj.addArcByCoordinates(xc,yc,v1(1),v1(2),v2(1),v2(2),extAngle>0);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addArcByCoordinates(xc,yc,v1(1),v1(2),v2(1),v2(2),extAngle>0);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = extendArcBySegment(obj, eIndex, extAngle, extAmplitude, seIndex)

            % set default start/end index as end
            if nargin < 5, seIndex = 1; end

            % get arc edge handle
            edgeHandle = obj.edges(eIndex).ptr;

            % check start/end point index: 0 when we extend inward, 1: for outward extension
            if seIndex == 0

                u = edgeHandle.getu1;
                if edgeHandle.direction
                    u = emdlab_g2d_rotatePoints(u,pi/2+extAngle);
                    v2 = edgeHandle.p1.getVector;
                    v1 = v2 - extAmplitude * u;
                else
                    u = emdlab_g2d_rotatePoints(u,-pi/2+extAngle);
                    v2 = edgeHandle.p1.getVector;
                    v1 = v2 - extAmplitude * u;
                end

            elseif seIndex == 1

                u = edgeHandle.getu2;
                if edgeHandle.direction
                    u = emdlab_g2d_rotatePoints(u,pi/2+extAngle);
                    v1 = edgeHandle.p2.getVector;
                    v2 = v1 + extAmplitude * u;
                else
                    u = emdlab_g2d_rotatePoints(u,-pi/2+extAngle);
                    v1 = edgeHandle.p2.getVector;
                    v2 = v1 + extAmplitude * u;
                end

            else

                error('Start/end index must be 0 or 1.');

            end

            if nargout == 0
                obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 1
                varargout{1} = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = extendArcByTangentSegment(obj, eIndex, seIndex, extAmplitude)

            edgeHandle = obj.edges(eIndex).ptr;

            % check start or index: 0 when we exten inward, 0: outward extension
            if seIndex == 0

                u = edgeHandle.getu1;
                if edgeHandle.direction
                    u = emdlab_g2d_rotatePoints(u,-pi/2);
                    v2 = edgeHandle.p1.getVector;
                    v1 = v2 + extAmplitude * u;
                else
                    u = emdlab_g2d_rotatePoints(u,pi/2);
                    v2 = edgeHandle.p1.getVector;
                    v1 = v2 + extAmplitude * u;
                end

            else

                u = edgeHandle.getu2;
                if edgeHandle.direction
                    u = emdlab_g2d_rotatePoints(u,pi/2);
                    v1 = edgeHandle.p2.getVector;
                    v2 = v1 + extAmplitude * u;
                else
                    u = emdlab_g2d_rotatePoints(u,-pi/2);
                    v1 = edgeHandle.p2.getVector;
                    v2 = v1 + extAmplitude * u;
                end

            end

            if nargout == 0
                obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 1
                varargout{1} = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = extendArcByPerpendicularSegment(obj, eIndex, seIndex, extAmplitude)

            % get arc handle
            edgeHandle = obj.edges(eIndex).ptr;

            % check start or index: 0 when we exten inward, 0: outward extension
            if seIndex == 0

                u = -edgeHandle.getu1;
                if edgeHandle.direction
                    u = emdlab_g2d_rotatePoints(u,-pi/2);
                    v2 = edgeHandle.p1.getVector;
                    v1 = v2 + extAmplitude * u;
                else
                    u = emdlab_g2d_rotatePoints(u,pi/2);
                    v2 = edgeHandle.p1.getVector;
                    v1 = v2 + extAmplitude * u;
                end

            else

                u = edgeHandle.getu2;
                if edgeHandle.direction
                    u = emdlab_g2d_rotatePoints(u,pi/2);
                    v1 = edgeHandle.p2.getVector;
                    v2 = v1 + extAmplitude * u;
                else
                    u = emdlab_g2d_rotatePoints(u,-pi/2);
                    v1 = edgeHandle.p2.getVector;
                    v2 = v1 + extAmplitude * u;
                end

            end

            if nargout == 0
                obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 1
                varargout{1} = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addSegmentByCoordinates(v1(1),v1(2),v2(1),v2(2));
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = extendArcByArc(obj, eIndex, xc, yc, extAngle, seIndex)

            % set default start/end index as end
            if nargin < 6, seIndex = 1; end

            % get edge handle
            edgeHandle = obj.edges(eIndex).ptr;

            % check start or end point of the segment for extension
            if seIndex == 0
                v1 = edgeHandle.p0.getVector;
            else
                v1 = edgeHandle.p1.getVector;
            end

            v2 = emdlab_g2d_rotatePoints(v1, extAngle, xc, yc);
            if nargout == 0
                obj.addArcByCoordinates(xc,yc,v1(1),v1(2),v2(1),v2(2),extAngle>0);
            elseif nargout == 1
                varargout{1} = obj.addArcByCoordinates(xc,yc,v1(1),v1(2),v2(1),v2(2),extAngle>0);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addArcByCoordinates(xc,yc,v1(1),v1(2),v2(1),v2(2),extAngle>0);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        %% edge edits
        function newEdgeIndex = splitEdge(obj, edgeIndex, splitRatio)

            if sum(splitRatio) > 1
                error('The summation of the split ratio mast be lower than one.');
            end
            eptr = obj.edges(edgeIndex).ptr;

            newEdgeIndex(1) = edgeIndex;

            switch class(eptr)
                case 'emdlab_g2d_segment'

                    p1Index = obj.getPointIndexByTag(eptr.p1.tag);
                    newp = eptr.p0.getVector;
                    vec = eptr.getUnitVector * eptr.getLength;
                    pnewIndex = [];
                    for i = 1:length(splitRatio)
                        newp = newp + splitRatio(i) * vec;
                        pnewIndex(end+1) = obj.addPoint(newp(1),newp(2));

                    end

                    eptr.p1 = obj.points(pnewIndex(1));
                    for i = 1:length(splitRatio)-1
                        newEdgeIndex(end+1) = obj.addSegment(pnewIndex(i),pnewIndex(i+1));
                    end

                    newEdgeIndex(end+1) = obj.addSegment(pnewIndex(end),p1Index);

            end

        end
        
        function removeEdge(obj, eIndex)

            % first remove all connected loops
            for lTag = obj.edges(eIndex).tags

            end

        end

        %% loop methods
        % adding a new loop to data base
        % this function returns loop index and loop handle
        function varargout = addLoop(obj, varargin)

            % get a loop class instance
            loopHandle = emdlab_g2d_loop;
            obj.loops(end+1) = loopHandle;

            % loop to assign edges and directions
            for i = 1:numel(varargin)

                for j = 1:numel(varargin{i})

                    if ischar(varargin{i}(j))

                        % get edge tag and remove spaces
                        edgeTag = strrep(varargin{i}(j), ' ', '');

                        if strcmpi(edgeTag(1), '-')
                            edgeIndex = obj.getEdgeIndexByTag(edgeTag(2:end));
                            loopHandle.addEdge(obj.edges(edgeIndex).ptr, false, -edgeIndex);
                        else
                            edgeIndex = obj.getEdgeIndexByTag(edgeTag);
                            loopHandle.addEdge(obj.edges(edgeIndex).ptr, true, edgeIndex);
                        end

                    else

                        edgeIndex = abs(varargin{i}(j));
                        loopHandle.addEdge(obj.edges(edgeIndex).ptr, varargin{i}(j)>0, varargin{i}(j));

                    end

                end

            end

            % generate loop tag
            loopIndex = numel(obj.loops);
            loopHandle.tag = ['l', num2str(loopIndex)];

            % add loop tag to connected edges
            for edgeIndex = abs(loopHandle.edgesIndexList)
                obj.edges(edgeIndex).tags(end+1) = loopHandle.tag;
            end

            if nargout == 1
                varargout{1} = loopIndex;
            elseif nargout == 2
                varargout{1} = loopIndex;
                varargout{2} = loopHandle;
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function removeLoop(obj, loopIndex)

            if ischar(loopIndex)
                loopIndex = obj.getLoopIndexByTag(loopIndex);
            end

            % first remove all connected faces to this loop
            for faceTag = obj.loops(loopIndex).tags
                obj.removeFace(faceTag);
            end

            obj.loops(loopIndex) = [];

        end

        function loopIndex = getLoopIndexByTag(obj, lTag)

            % check for existance of already defined loop in data base
            for i = 1:numel(obj.loops)

                if strcmp(obj.loops(i).tag,lTag)

                    loopIndex = i;
                    return;

                end

            end

            error('Loop was not found.');

        end

        %% face methods
        % adding a new face to data base
        % this function returns the face handle
        function varargout = addFace(obj, faceName, varargin)

            % get face class instance
            faceHandle = emdlab_g2d_face;
            faceHandle.tag = faceName;
            faceHandle.color = rand(1,3);
            obj.faces(end+1) = faceHandle;

            for i = 1:numel(varargin)

                for j = 1:numel(varargin{i})
                    faceHandle.addLoop(obj.loops(varargin{i}(j)));
                    obj.loops(varargin{i}(j)).tags(end+1) = faceName;
                end

            end            

            if nargout == 1
                varargout{1} = faceHandle;
            elseif nargout > 1
                error('The number of output arguments is too hight');
            end

        end

        function removeFace(obj, faceIndex)

            if ischar(faceIndex)
                faceIndex = obj.getFaceIndexByTag(faceIndex);
            end
            obj.faces(faceIndex) = [];
            
        end

        function faceIndex = getFaceIndexByTag(obj, fTag)

            % check for existance of already defined face in data base
            for i = 1:numel(obj.faces)

                if strcmp(obj.faces(i).tag,fTag)

                    faceIndex = i;
                    return;

                end

            end

            error('Face was not found.');

        end

        function setFaceColor(obj, faceTag, R, G, B)
            obj.faces(obj.getFaceIndexByTag(faceTag)).color = [R,G,B]/255;
        end

        %% mesh generation methods
        % generate triangular mesh for geometry
        function m = generateMesh(obj, meshGenerator)

            % default mesh generator
            if nargin<2
                meshGenerator = 'mm';
            end

            if strcmpi(meshGenerator, 'gmsh')
                obj.write_geo_file;
                m = obj.read_msh_file;
                return;
            end

            % get an instance of mesh data base
            m = emdlab_m2d_tmdb;

            % add mesh zones
            for i = 1:numel(obj.faces)
                m.addMeshZone(obj.faces(i).tag, obj.faces(i).getMesh(meshGenerator));
                m.mzs.(obj.faces(i).tag).color = obj.faces(i).color;
            end

        end

        %% visualization methos
        % show the geometry sketch
        function varargout = showSketch(obj, showTags, showWFM)

            if nargin<2
                showTags = true;
                showWFM = false;
            elseif nargin<3
                showWFM = false;
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
                plot(p(:,1), p(:,2), 's', 'LineWidth', 1.5, 'MarkerEdgeColor','k');
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
                if showTags
                    
                    edgeTags = cell(1,Ne);
                    for i = 1:Ne
                        edgeTags{i} = obj.edges(i).tag;
                    end
                    text(c(:,1), c(:,2), edgeTags, 'BackgroundColor', 'w', ...
                        'HorizontalAlignment','center','VerticalAlignment','middle');
                end
                if showWFM
                    plot(v(:,1), v(:,2), 'o', 'color', 'k', 'markersize',5, 'markerfacecolor','k');
                end
            end

            set(gca, 'clipping', 'off');

            axis off equal
            zoom on
            grid on
            drawnow;

            if nargout == 1
                varargout{1} = f;
            end

        end

        function showFaces(obj)

            m = obj.generateMesh('mm');
            m.showg;

        end
        %% adding primitive loops
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
                obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addCornerRectangleLoop(obj, x0, y0, w, h)

            p1Index = obj.addPoint(x0,y0);
            p2Index = obj.addPoint(x0+w,y0);
            p3Index = obj.addPoint(x0+w,y0+h);
            p4Index = obj.addPoint(x0,y0+h);

            e1Index = obj.addSegment(p1Index, p2Index);
            e2Index = obj.addSegment(p2Index, p3Index);
            e3Index = obj.addSegment(p3Index, p4Index);
            e4Index = obj.addSegment(p4Index, p1Index);

            if nargout == 0
                obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addCenterRectangleLoop(obj, x0, y0, w, h)

            p1Index = obj.addPoint(x0-w/2,y0-h/2);
            p2Index = obj.addPoint(x0+w/2,y0-h/2);
            p3Index = obj.addPoint(x0+w/2,y0+h/2);
            p4Index = obj.addPoint(x0-w/2,y0+h/2);

            e1Index = obj.addSegment(p1Index, p2Index);
            e2Index = obj.addSegment(p2Index, p3Index);
            e3Index = obj.addSegment(p3Index, p4Index);
            e4Index = obj.addSegment(p4Index, p1Index);

            if nargout == 0
                obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = add3PointCornerRectangleLoop(obj, x0, y0, x1, y1, x2, y2)

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

        function varargout = add3PointCenterRectangleLoop(obj, x0, y0, w, h)

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

        function varargout = addParallelogramLoop(obj, x0, y0, w, h)

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
                obj.addLoop(e1Index, e2Indexe);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, e2Indexe);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, e2Indexe);
            elseif nargout > 2
                error('The number of output arguments is too high.');
            end

        end

        function varargout = addAnnularSectorLoop(obj, Ri, Ro, Theta_1, Theta_2, xc, yc)

            % set default center
            if nargin < 6
                xc = 0;
                yc = 0;
            end

            oIndex = obj.addPoint(xc, yc);
            p1Index = obj.addPoint(xc + Ri * cos(Theta_1), yc + Ri * sin(Theta_1));
            p2Index = obj.addPoint(xc + Ro * cos(Theta_1), yc + Ro * sin(Theta_1));
            p3Index = obj.addPoint(xc + Ro * cos(Theta_2), yc + Ro * sin(Theta_2));
            p4Index = obj.addPoint(xc + Ri * cos(Theta_2), yc + Ri * sin(Theta_2));

            e1Index = obj.addSegment(p1Index, p2Index);
            e2Index = obj.addArc(oIndex, p2Index, p3Index, 1);
            e3Index = obj.addSegment(p3Index, p4Index);
            e4Index = obj.addArc(oIndex, p4Index, p1Index, 0);

            if nargout == 0
                obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 1
                varargout{1} = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
            elseif nargout == 2
                [varargout{1},varargout{2}] = obj.addLoop(e1Index, e2Index, e3Index, e4Index);
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
            e_indices = cell(1,Nx);
            for i = 1:Nx
                e_indices{i} = obj.addSegment(p_indices(i),p_indices(i+1));
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

            % set mesh size at points for Gmsh
            for i = 1:numel(obj.points)
                obj.points(i).meshSize = mLength;
            end

        end

        function setLoopMeshMaxLength(obj, loopIndex, mLength)

            % get loop handle
            loopHandle = obj.loops(loopIndex);

            for i = 1:numel(loopHandle.edges)
                eptr = loopHandle.edges{i};
                eptr.setMaxLength(mLength);
                switch class(eptr)
                    case 'emdlab_g2d_segment'
                        eptr.p0.meshSize = mLength;
                        eptr.p1.meshSize = mLength;
                    case 'emdlab_g2d_arc'
                        eptr.p1.meshSize = mLength;
                        eptr.p2.meshSize = mLength;
                end
            end

        end

        function setMeshLengthByRadialFunction(obj, fHandle)

            for i = 1:numel(obj.edges)

                if isa(obj.edges(i).ptr, 'emdlab_g2d_segment')
                    obj.edges(i).ptr.setL1L2(fHandle(obj.edges(i).ptr.p0.norm()), fHandle(obj.edges(i).ptr.p1.norm()))
                elseif  isa(obj.edges(i).ptr, 'emdlab_g2d_arc')
                    obj.edges(i).ptr.setMaxLength(fHandle(obj.edges(i).ptr.p1.norm()));
                elseif  isa(obj.edges(i).ptr, 'emdlab_g2d_spline')
                    obj.edges(i).ptr.setL1L2(fHandle(obj.edges(i).ptr.pts(1).norm()), fHandle(obj.edges(i).ptr.pts(end).norm()));
                end

            end

            % set mesh size at points for Gmsh
            for i = 1:numel(obj.points)
                obj.points(i).meshSize = fHandle(obj.points(i).norm());
                if isnan(obj.points(i).meshSize), obj.points(i).meshSize = 1; end
            end

        end

        function setMeshLengthByYFunction(obj, fHandle)

            for i = 1:numel(obj.edges)

                if isa(obj.edges(i).ptr, 'emdlab_g2d_segment')
                    obj.edges(i).ptr.setL1L2(fHandle(obj.edges(i).ptr.p0.y), fHandle(obj.edges(i).ptr.p1.y))
                elseif  isa(obj.edges(i).ptr, 'emdlab_g2d_arc')
                    obj.edges(i).ptr.setMaxLength(fHandle(obj.edges(i).ptr.p1.y));
                elseif  isa(obj.edges(i).ptr, 'emdlab_g2d_spline')
                    obj.edges(i).ptr.setL1L2(fHandle(obj.edges(i).ptr.pts(1).y), fHandle(obj.edges(i).ptr.pts(end).y));
                end

            end

            % set mesh size at points for Gmsh
            for i = 1:numel(obj.points)
                obj.points(i).meshSize = fHandle(obj.points(i).y);
                if isnan(obj.points(i).meshSize), obj.points(i).meshSize = 1; end
            end

        end

        % interact with Gmsh software
        function write_geo_file(obj)

            % define a new geo file
            fid = fopen("C:\emdlab-win64\tmp\emdlab_gmsh_geoFile.geo", 'w');
            fprintf(fid, 'SetFactory("OpenCASCADE");\n');

            % add points
            for i = 1:numel(obj.points)
                fprintf(fid, 'Point(%d) = {%.16f, %.16f, 0, %f};\n', i, obj.points(i).x, obj.points(i).y, obj.points(i).meshSize);
            end

            % add edges
            for i = 1:numel(obj.edges)

                if isa(obj.edges(i).ptr, 'emdlab_g2d_segment')

                    fprintf(fid, 'Line(%d) = {%d, %d};\n', i, obj.getPointIndexByTag(obj.edges(i).ptr.p0.tag), ...
                        obj.getPointIndexByTag(obj.edges(i).ptr.p1.tag));

                elseif isa(obj.edges(i).ptr, 'emdlab_g2d_arc')

                    fprintf(fid, 'Circle(%d) = {%d, %d, %d};\n', i, obj.getPointIndexByTag(obj.edges(i).ptr.p1.tag), ...
                        obj.getPointIndexByTag(obj.edges(i).ptr.p0.tag), obj.getPointIndexByTag(obj.edges(i).ptr.p2.tag));

                elseif isa(obj.edges(i).ptr, 'emdlab_g2d_spline')

                    pointsList = zeros(1,numel(obj.edges(i).ptr.pts));
                    for j = 1:numel(obj.edges(i).ptr.pts)
                        pointsList(j) = obj.getPointIndexByTag(obj.edges(i).ptr.pts(j).tag);
                    end
                    pointsList = join(string(pointsList), ", ");
                    fprintf(fid, 'Spline(%d) = {%s};\n', i, pointsList);

                end

            end

            % add loops
            for i = 1:numel(obj.loops)

                fprintf(fid, 'Curve Loop(%d) = {%s};\n',i, join(string(obj.loops(i).edgesIndexList), ', '));

            end

            % add faces
            for i = 1:numel(obj.faces)

                tmp_str = num2str(obj.getLoopIndexByTag(obj.faces(i).loops(1).tag));
                for j = 2:length(obj.faces(i).loops)
                    tmp_str = [tmp_str , ', ' , num2str(obj.getLoopIndexByTag(obj.faces(i).loops(j).tag))];
                end

                fprintf(fid, 'Plane Surface(%d) = {%s};\n', i, tmp_str);
                fprintf(fid, 'Physical Surface(%d) = {%d};\n', i, i);

            end

            fclose(fid);

        end

        function m = read_msh_file(obj)

            % run gmsh via matlab            
            pyCodePath = "C:\\emdlab-win64\\py-files\\gmsh\\emdlab_gmsh_runGeoSaveMsh2D.py";

            [~,~] = system(char('"' + obj.pyPath + '"' + " " + '"' + pyCodePath+ '"'));

            % read generated mesh;
            emdlab_gmsh_mshFile;

            % get an instance of mesh data base
            m = emdlab_m2d_tmdb;

            nodes = msh.POS(:,1:2);
            Np = size(nodes,1);

            % add faces
            for i = 1:numel(obj.faces)

                cl = msh.TRIANGLES(msh.TRIANGLES(:,4) == i, 1:3);
                index = unique(cl(:));
                index = sort(index);
                xpoints = nodes(index,:);
                pindex = zeros(Np,1);
                pindex(index) = 1:size(xpoints,1);
                cl = pindex(cl);

                p21 = xpoints(cl(:,2),:) - xpoints(cl(:,1),:);
                p31 = xpoints(cl(:,3),:) - xpoints(cl(:,1),:);
                index = (p21(:,1).*p31(:,2) - p21(:,2).*p31(:,1)) < 0;
                cl(index,:) = cl(index,[1,3,2]);

                m.addMeshZone(obj.faces(i).tag, emdlab_m2d_tmz(cl, xpoints));
                m.mzs.(obj.faces(i).tag).color = obj.faces(i).color;

            end

        end

        function extrudeAndSaveStepSTL(obj, faceName, z1, z2)

            obj.write_geo_file;
            index = obj.getFaceIndexByTag(faceName);

            % define a new geo file
            fid = fopen("C:\emdlab-win64\tmp\emdlab_gmsh_geoFile.geo", 'a');
            fprintf(fid, "Extrude {0, 0, %.16f} {Surface{%s};}\n", z2-z1, join(string(1:numel(obj.faces)),','));
            fprintf(fid, "Translate {0, 0, %.16f} {Volume{%d};}\n", z1, index);
            fprintf(fid, "Recursive Delete { Volume{%s};}\n", join(string(setdiff(1:numel(obj.faces), index)),','));
            fprintf(fid, 'Coherence;\n');
            fclose(fid);

            % run gmsh via matlab            
            pyCodePath = "C:\\emdlab-win64\\py-files\\emdlab_gmsh_runGeoSaveStep.py";

            [~,~] = system(char('"' + obj.pyPath + '"' + " " + '"' + pyCodePath+ '"'));

            stpPath = "C:\emdlab-win64\geometry\step\emdlab_g3d_stepFile.step";
            copyfile(stpPath, cd + "\" + string(faceName) + ".step")


        end

    end

end