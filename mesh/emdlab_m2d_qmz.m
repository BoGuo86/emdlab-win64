% developer: https://ComProgExpert.com, Ali Jamali-Fard
% 2D quadrilateral mesh zone

classdef emdlab_m2d_qmz <  handle & emdlab_g2d_constants & matlab.mixin.Copyable
    
    properties (SetAccess = private)

        % mesh nodes
        nodes (:,2) double;

        % mesh connectivity list
        cl (:,4) double;

        % mesh elements
        elements (:,4) double;

        % unique edges
        edges (:,2) double;

        % list of boundary edges
        bedges (:,1) double;
        
        acoefs (1,1) double;
        bcoefs (1,1) double;

    end
    properties

        % zone index
        zi (1,1) double;

        % local to global node index
        ni_l2g (:,1) double;
        ei_l2g (:,1) double;

        % material of zone
        material char = 'air';

        % mesh zone color
        color = 'c';

        % mesh zone properties: differs in differents solvers
        props (1,1) struct;

    end

    properties (Access = private)

        % element area
        ea (:,1) double;

        % mesh zone area
        area (1,1) double;

        % states
        isDataSetted (1,1) logical = false;
        is_ea_Evaluated (1,1) logical = false;
        is_area_Evaluated (1,1) logical = false;

    end

    properties (Dependent = true)
        
        % number of mesh zone nodes
        Nn(1, 1) double;

        % number of mesh zone elements
        Ne(1, 1) double;
        
    end

    methods
        %% Constructor and Destructor
        function obj = emdlab_m2d_qmz(cl, nodes)
            obj.nodes = nodes;
            obj.cl = cl;
            obj = obj.setdata;
        end

        function y = get.Nn(obj)
            y = size(obj.nodes, 1);
        end

        function y = get.Ne(obj)
            y = size(obj.cl, 1);
        end
        
        %% Topological Functions
        % setting needed data
        function obj = setdata(obj)

            % check if already data is set
            if obj.isDataSetted, return; end

            % first edge of each quadrilateral
            e1 = obj.cl(:,[1,2]);
            % second edge of each quadrilateral
            e2 = obj.cl(:,[2,3]);
            % third edge of each quadrilateral
            e3 = obj.cl(:,[3,4]);
            % forth edge of each quadrilateral
            e4 = obj.cl(:,[4,1]);

            % sorting for lower index
            [e1,s1] = sort(e1, 2);
            [e2,s2] = sort(e2, 2);
            [e3,s3] = sort(e3, 2);
            [e4,s4] = sort(e4, 2);

            % specefying changed edge index
            s1 = s1(:,1) == 2;
            s2 = s2(:,1) == 2;
            s3 = s3(:,1) == 2;
            s4 = s4(:,1) == 2;

            % unification of edges
            [obj.edges, ~, ic] = unique([e1; e2; e3; e4],'rows');

            % getting number of elements
            ne = obj.Ne;

            % getting index of edge corresponding to each elements
            e1 = ic(1:ne);
            e2 = ic(1+ne:2*ne);
            e3 = ic(1+2*ne:3*ne);
            e4 = ic(1+3*ne:4*ne);

            % specefying boundary edges
            obj.bedges = sparse([e1,e2,e3,e4],ones(4*ne,1),ones(4*ne,1));
            obj.bedges = full(obj.bedges==1);

            % specefying trace direction
            e1(s1) = -e1(s1);
            e2(s2) = -e2(s2);
            e3(s3) = -e3(s3);
            e4(s4) = -e4(s4);

            % element matrix
            obj.elements = [e1,e2,e3,e4];

            % change states
            obj.isDataSetted = true;

        end
        
        %% Mesh Visiualization
        function varargout = showm(obj)
            
            f = emdlab_r2d_mesh();
            ax = axes(f);
            f.Name = ['[Mesh Zone][', 'Nn = ', num2str(obj.Nn), '][Ne = ', num2str(obj.Ne), ']'];

            patch('Faces',obj.cl(:,1:4),'Vertices',obj.nodes,'FaceColor',...
                obj.color,'EdgeColor','k','parent',ax);
            
            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout == 2
                varargout{1} = f;
                varargout{2} = ax;
            elseif nargout > 2
                error('Too many output argument.');
            end

        end
        function varargout = showwf(obj)


            f = emdlab_r2d_mesh();
            ax = axes(f);

            patch('Faces',obj.edges(find(obj.bedges),:),'Vertices',obj.nodes,...
                'FaceColor','none','EdgeColor','k', 'parent', ax);
            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout == 2
                varargout{1} = f;
                varargout{2} = ax;
            elseif nargout > 2
                error('Too many output argument.');
            end
        end
        
        %% Tools Functions
        function obj = moveNodes(obj,MovTol)
            if nargin<2
                MovTol = 1e-3;
            end
            % connectivity matrix for nodes
            Con = sparse(obj.edges(:,1),obj.edges(:,2),...
                ones(size(obj.edges,1),1),obj.Nn,obj.Nn);
            Con = Con + Con';
            % loop for movments
            inodes = obj.getinodes;
            % weight matrix
            weight = diag(1./sum(Con(inodes,:),2));
            for iter = 1:100
                % getting position of new nodes
                pnew = Con(inodes,:)*obj.nodes;
                pnew = weight*pnew;
                % evaluation of movments
                Mov = sqrt(sum((obj.nodes(inodes,:)-pnew).^2,2));
                disp(sum(Mov));
                obj.nodes(inodes,:) = pnew;
                % check for movment tolerance
                if Mov < MovTol
                    disp(iter);
                    break
                end
            end
        end
        function y = getbnodes(obj)
            % getting index of boundary nodes
            y = obj.edges(obj.bedges,:);
            y = unique(y(:));
        end
        function y = getinodes(obj)
            % getting index of inner nodes
            y = obj.getbnodes;
            y = setdiff((1:obj.Nn)',y);
        end
        function strefine(obj)
            % number of nodes and elements in old mesh
            NnOld = obj.Nn;
            NeOld = obj.Ne;
            NedOld = size(obj.edges, 1);
            % nodes of new mesh
            obj.nodes = [obj.nodes;...
                (obj.nodes(obj.edges(:,1),:)+obj.nodes(obj.edges(:,2),:))/2;
                (obj.nodes(obj.cl(:,1),:)+obj.nodes(obj.cl(:,2),:)+obj.nodes(obj.cl(:,3),:)+obj.nodes(obj.cl(:,4),:))/4];
            % index of nodes on old edges
            index = [abs(obj.elements),(1:NeOld)'];
            % new connctivity list
            obj.cl = [obj.cl(:,1),index(:,1)+NnOld,index(:,5)+NnOld+NedOld,index(:,4)+NnOld
                obj.cl(:,2),index(:,2)+NnOld,index(:,5)+NnOld+NedOld,index(:,1)+NnOld
                obj.cl(:,3),index(:,3)+NnOld,index(:,5)+NnOld+NedOld,index(:,2)+NnOld
                obj.cl(:,4),index(:,4)+NnOld,index(:,5)+NnOld+NedOld,index(:,3)+NnOld];
            % setting data of new mesh
            obj.setdata;
        end
    
    end
end
