classdef emdlab_m3d_ttmdb < handle & emdlab_g2d_constants
    
    properties (SetAccess = private)
        
        % mesh zones
        mzs
        % Mesh Nodes
        nodes(:, 3) double;
        % Mesh Connectivity List
        cl(:, 4) double;
        % Unique Mesh Edges
        edges
        % Unique Mesh Facets
        facets
        % Boundary Facets
        bfacets
        % Tetrahedral Elements
        elements(:, 5) double;
        % stiffness, force, mass and connectivity matrix
        Ke
        Fe
        Me
        % jacobian inverse transpose
        JIT
        % element zone index
        ezi
        % global elements volume
        gev
        % gradient
        gphix
        gphiy
        gphiz
        % materials
        mts
        % auxiliary stored matricies
        mtcs
        % named selections
        facetNamedSelections
        
    end
    
    properties (Dependent = true)
        
        % Number of nodes
        Nn(1,1) double;
        % Number of elements
        Ne(1,1) double;
        % number of mesh zones
        Nmzs(1,1) double;
        
    end
    
    properties (Access = private)
        
        % element type
        etype
        % states
        isGlobalMeshGenerated = false;
        is_JIT_Evaluated = false;
        isKeFe_TTL4_Evaluated = false;
        isKexyz9Fe_Evaluated = false;
        isKexyz1_TTL4_Evaluated = false;
        isKexyz3_TTL4_Evaluated = false;
        
    end
    
    methods
        %% Initialization
        function obj = emdlab_m3d_ttmdb(varargin)
            
            % default values
            obj.mtcs = struct();
            obj.mts = struct();
            obj.mzs = struct();
            
        end
        
        %% FEM Preparation
        function ggmesh(obj)
            if obj.isGlobalMeshGenerated, return; end
            % generation of initial mesh
            Nn_tmp = 0;
            Ne_tmp = 0;
            mzNames = fieldnames(obj.mzs);
            
            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                Nn_tmp = Nn_tmp + mzptr.Nn;
                Ne_tmp = Ne_tmp + mzptr.Ne;
            end
            
            % initialization of nodes and elements
            obj.nodes = zeros(Nn_tmp, 3);
            obj.cl = zeros(Ne_tmp, 4);
            obj.elements = zeros(Ne_tmp, 5);
            nindex = 0;
            eindex = 0;
            
            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                % insertion of nodes
                obj.nodes(1 + nindex:mzptr.Nn + nindex, :) = mzptr.nodes;
                % insertion of elements
                obj.cl(1 + eindex:mzptr.Ne + eindex, :) = mzptr.cl + nindex;
                % specefying zone index
                obj.elements(1 + eindex:mzptr.Ne + eindex, 5) = i;
                mzptr.zi = i;
                nindex = nindex + mzptr.Nn;
                eindex = eindex + mzptr.Ne;
            end
            
            [obj.nodes, ~, ic] = uniquetol(obj.nodes, obj.gleps, 'ByRows', true);
            obj.cl = ic(obj.cl);
            % setting l2g
            nindex = 0;
            
            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                mzptr.l2g = ic(nindex + 1:nindex + mzptr.Nn);
                nindex = nindex + mzptr.Nn;
            end
            
            obj.setdata;
            obj.evalezi;
            obj.isGlobalMeshGenerated = true;
        end
        
        function evalezi(obj)
            obj.ezi = zeros(obj.Ne, obj.Nmzs, 'logical');
            
            for i = 1:obj.Nmzs
                obj.ezi(:, i) = obj.elements(:, 5) == i;
            end
            
        end
        
        function evalJIT(obj)
            % check states
            if obj.is_JIT_Evaluated, return; end
            tic, disp('-------------------------------------------------------');
            % prerequisite
            obj.ggmesh;
            % gev: global element volume
            obj.gev = zeros(1, obj.Ne);
            mzNames = fieldnames(obj.mzs);
            
            for i = 1:obj.Nmzs
                obj.gev(obj.ezi(:, obj.mzs.(mzNames{i}).zi)) = ...
                    obj.mzs.(mzNames{i}).getVolumeOfElements;
            end
            
            % evaluation of jacobian inverse using d1 element data
            obj.JIT = ttmdbc_evalJIT(obj.cl, obj.nodes, obj.gev);
            % change states
            obj.is_JIT_Evaluated = true;
            disp('Evaluation of JIT completed.');
            toc, disp('-------------------------------------------------------');
        end
        
        function evalKeFe_TTL4(obj)
            if obj.isKeFe_TTL4_Evaluated, return; end
            obj.ggmesh;
            obj.evalJIT;
            % getting number of elements
            xNe = obj.Ne;
            % getting data of reference element
            tic, disp('-------------------------------------------------------');
            % evaluation of grad phi i
            [obj.gphix, obj.gphiy, obj.gphiz] = ttmdbc_evalGphixyz_TTL4(obj.JIT);
            % evaluation of Ke
            obj.Ke = zeros(10, xNe);
            temp = 0;
            
            for i = 1:4
                
                for j = 1:i
                    temp = temp + 1;
                    obj.Ke(temp, :) = ...
                        obj.gphix(i, :) .* obj.gphix(j, :) + ...
                        obj.gphiy(i, :) .* obj.gphiy(j, :) + ...
                        obj.gphiz(i, :) .* obj.gphiz(j, :);
                end
                
            end
            
            % multiplying by triangle areas
            obj.Ke = obj.Ke * sparse(1:xNe, 1:xNe, obj.gev);
            % elemental force matrix
            obj.Fe = repmat((obj.gev / 4), 4, 1);
            disp('Evaluation of [Ke] and [Fe] completed.');
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKeFe_TTL4_Evaluated = true;
        end
        
        function evalKeFe(obj, etype)
            xNe = obj.Ne;
            obj.etype = etype;
            % getting data of reference element
            tic
            % evaluation of jacobian inverse using d1 element data
            edata = getedata('TTL4');
            % x, y and z coordinate of points
            xp = obj.nodes(:, 1);
            yp = obj.nodes(:, 2);
            zp = obj.nodes(:, 3);
            % point coordinate of each triangle nodes
            xp = xp(obj.cl');
            yp = yp(obj.cl');
            zp = zp(obj.cl');
            % evaluation of a, b and c coefficient
            acoefs = edata.M \ xp;
            bcoefs = edata.M \ yp;
            ccoefs = edata.M \ zp;
            % setting elements volume
            obj.gev = zeros(1, obj.Ne);
            mznames = fieldnames(obj.mzs);
            
            for i = 1:obj.Nmzs
                obj.gev(obj.ezi(:, obj.mzs.(mznames{i}).zi)) = ...
                    obj.mzs.(mznames{i}).ev;
            end
            
            % Jacobiab transpose matrix of each elements
            obj.JT = [acoefs(2, :); acoefs(3, :); acoefs(4, :); ...
                bcoefs(2, :); bcoefs(3, :); bcoefs(4, :); ...
                ccoefs(2, :); ccoefs(3, :); ccoefs(4, :); ];
            [i, j] = getij(3, xNe);
            obj.JT = sparse(i, j, obj.JT(:));
            % getting specefied edata
            edata = getedata(obj.etype);
            
            switch obj.etype
                case 'TTL4'
                    gphi = obj.JT \ repmat(edata.GG, xNe, 1);
                    % inserting in gradx and grady phi i
                    obj.gphix = (gphi(1:3:end, :))';
                    obj.gphiy = (gphi(2:3:end, :))';
                    obj.gphiz = (gphi(3:3:end, :))';
                    % evaluation of Ke
                    obj.Ke = zeros(10, xNe);
                    temp = 0;
                    
                    for i = 1:4
                        
                        for j = 1:i
                            temp = temp + 1;
                            obj.Ke(temp, :) = ...
                                obj.gphix(i, :) .* obj.gphix(j, :) + ...
                                obj.gphiy(i, :) .* obj.gphiy(j, :) + ...
                                obj.gphiz(i, :) .* obj.gphiz(j, :);
                        end
                        
                    end
                    
                    % multiplying by triangle areas
                    obj.Ke = obj.Ke * sparse(1:xNe, 1:xNe, obj.gev);
                    % elemental force matrix
                    obj.Fe = repmat((obj.gev / 4), 4, 1);
                otherwise
                    error('Element type does not defined.');
            end
            
            disp('************************************')
            disp('Calculation of Ke, Fe and C completed.')
            toc
        end
        
        function evalKeFeTTL4_Nodal(obj)
            obj.evalJIT;
            tic, disp('-------------------------------------------------------');
            xNe = obj.Ne;
            % getting data of reference element
            tic
            % evaluation of grad phi i
            [obj.gphix, obj.gphiy, obj.gphiz] = ttmdbc_evalGphixyz_TTL4(obj.JIT);
            % evaluation of Ke
            obj.Ke.K11 = zeros(10, xNe);
            obj.Ke.K22 = zeros(10, xNe);
            obj.Ke.K33 = zeros(10, xNe);
            obj.Ke.K21 = zeros(16, xNe);
            obj.Ke.K31 = zeros(16, xNe);
            obj.Ke.K32 = zeros(16, xNe);
            temp = 0;
            
            for i = 1:4
                
                for j = 1:i
                    temp = temp + 1;
                    obj.Ke.K11(temp, :) = obj.gphiy(i, :) .* obj.gphiy(j, :) + obj.gphiz(i, :) .* obj.gphiz(j, :);
                    obj.Ke.K22(temp, :) = obj.gphix(i, :) .* obj.gphix(j, :) + obj.gphiz(i, :) .* obj.gphiz(j, :);
                    obj.Ke.K33(temp, :) = obj.gphiy(i, :) .* obj.gphiy(j, :) + obj.gphix(i, :) .* obj.gphix(j, :);
                end
                
            end
            
            temp = 0;
            
            for i = 1:4
                
                for j = 1:4
                    temp = temp + 1;
                    obj.Ke.K21(temp, :) = -1 * (obj.gphix(i, :) .* obj.gphiy(j, :));
                    obj.Ke.K31(temp, :) = -1 * (obj.gphix(i, :) .* obj.gphiz(j, :));
                    obj.Ke.K32(temp, :) = -1 * (obj.gphiy(i, :) .* obj.gphiz(j, :));
                end
                
            end
            
            % multiplying by triangle areas
            obj.Ke.K11 = obj.Ke.K11 * sparse(1:xNe, 1:xNe, obj.gev);
            obj.Ke.K22 = obj.Ke.K22 * sparse(1:xNe, 1:xNe, obj.gev);
            obj.Ke.K33 = obj.Ke.K33 * sparse(1:xNe, 1:xNe, obj.gev);
            obj.Ke.K21 = obj.Ke.K21 * sparse(1:xNe, 1:xNe, obj.gev);
            obj.Ke.K31 = obj.Ke.K31 * sparse(1:xNe, 1:xNe, obj.gev);
            obj.Ke.K32 = obj.Ke.K32 * sparse(1:xNe, 1:xNe, obj.gev);
            % elemental force matrix
            obj.Fe = repmat((obj.gev / 4), 4, 1);
            disp('Evaluation of [Ke] and [Fe] completed.')
            toc, disp('-------------------------------------------------------');
        end
        
        function evalKexyz1_TTL4(obj)
            if obj.isKexyz1_TTL4_Evaluated, return; end
            obj.ggmesh;
            obj.evalJIT;
            tic, disp('-------------------------------------------------------');
            % evaluation of grad phi i
            [obj.gphix, obj.gphiy, obj.gphiz] = ttmdbc_evalGphixyz_TTL4(obj.JIT);
            % evaluation of Ke
            obj.Ke = zeros(10, obj.Ne);
            temp = 0;
            
            for i = 1:4
                
                for j = 1:i
                    temp = temp + 1;
                    obj.Ke(temp, :) = (obj.gphix(i, :) .* obj.gphix(j, :) + ...
                        obj.gphiy(i, :) .* obj.gphiy(j, :) + ...
                        obj.gphiz(i, :) .* obj.gphiz(j, :)) .* obj.gev;
                end
                
            end
            
            disp('Evaluation of [Ke] and [Fe] completed.')
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKexyz1_TTL4_Evaluated = true;
        end
        
        function evalKexyz3_TTL4(obj)
            if obj.isKexyz3_TTL4_Evaluated, return; end
            obj.evalJIT;
            tic, disp('-------------------------------------------------------');
            % evaluation of grad phi i
            [obj.gphix, obj.gphiy, obj.gphiz] = ttmdbc_evalGphixyz_TTL4(obj.JIT);
            % evaluation of Ke
            obj.Ke.xx = zeros(10, obj.Ne);
            obj.Ke.yy = zeros(10, obj.Ne);
            obj.Ke.zz = zeros(10, obj.Ne);
            temp = 0;
            
            for i = 1:4
                
                for j = 1:i
                    temp = temp + 1;
                    obj.Ke.xx(temp, :) = obj.gphix(i, :) .* obj.gphix(j, :) .* obj.gev;
                    obj.Ke.yy(temp, :) = obj.gphiy(i, :) .* obj.gphiy(j, :) .* obj.gev;
                    obj.Ke.zz(temp, :) = obj.gphiz(i, :) .* obj.gphiz(j, :) .* obj.gev;
                end
                
            end
            
            disp('Evaluation of [Ke] and [Fe] completed.')
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKexyz3_TTL4_Evaluated = true;
        end
        
        function evalKexyz9FeTTL4(obj)
            if obj.isKexyz9Fe_Evaluated, return; end
            obj.evalJIT;
            tic, disp('-------------------------------------------------------');
            xNe = obj.Ne;
            % getting data of reference element
            tic
            % evaluation of grad phi i
            [obj.gphix, obj.gphiy, obj.gphiz] = ttmdbc_evalGphixyz_TTL4(obj.JIT);
            % evaluation of Ke
            obj.Ke.Kxx = zeros(10, xNe);
            obj.Ke.Kyy = zeros(10, xNe);
            obj.Ke.Kzz = zeros(10, xNe);
            obj.Ke.Kyx = zeros(16, xNe);
            obj.Ke.Kzx = zeros(16, xNe);
            obj.Ke.Kzy = zeros(16, xNe);
            temp = 0;
            
            for i = 1:4
                
                for j = 1:i
                    temp = temp + 1;
                    obj.Ke.Kxx(temp, :) = obj.gphix(i, :) .* obj.gphix(j, :) .* obj.gev;
                    obj.Ke.Kyy(temp, :) = obj.gphiy(i, :) .* obj.gphiy(j, :) .* obj.gev;
                    obj.Ke.Kzz(temp, :) = obj.gphiz(i, :) .* obj.gphiz(j, :) .* obj.gev;
                end
                
            end
            
            temp = 0;
            
            for i = 1:4
                
                for j = 1:4
                    temp = temp + 1;
                    obj.Ke.Kyx(temp, :) = obj.gphix(i, :) .* obj.gphiy(j, :) .* obj.gev;
                    obj.Ke.Kzx(temp, :) = obj.gphix(i, :) .* obj.gphiz(j, :) .* obj.gev;
                    obj.Ke.Kzy(temp, :) = obj.gphiy(i, :) .* obj.gphiz(j, :) .* obj.gev;
                end
                
            end
            
            % elemental force matrix
            obj.Fe = repmat((obj.gev / 4), 4, 1);
            disp('Evaluation of [Ke] and [Fe] completed.')
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKexyz9Fe_Evaluated = true;
        end
        
        function obj = evalKe6Fe(obj, etype)
            xNe = obj.Ne;
            obj.etype = etype;
            % getting data of reference element
            tic
            % evaluation of jacobian inverse using d1 element data
            edata = getedata('TTL4');
            % x, y and z coordinate of points
            xp = obj.nodes(:, 1);
            yp = obj.nodes(:, 2);
            zp = obj.nodes(:, 3);
            % point coordinate of each triangle nodes
            xp = xp(obj.cl');
            yp = yp(obj.cl');
            zp = zp(obj.cl');
            % evaluation of a, b and c coefficient
            acoefs = edata.M \ xp;
            bcoefs = edata.M \ yp;
            ccoefs = edata.M \ zp;
            % setting elements volume
            obj.gev = zeros(1, obj.Ne);
            mznames = fieldnames(obj.mzs);
            
            for i = 1:obj.Nmzs
                obj.gev(obj.ezi(:, obj.mzs.(mznames{i}).zi)) = ...
                    obj.mzs.(mznames{i}).ev;
            end
            
            % Jacobiab transpose matrix of each elements
            obj.JT = [acoefs(2, :); acoefs(3, :); acoefs(4, :); ...
                bcoefs(2, :); bcoefs(3, :); bcoefs(4, :); ...
                ccoefs(2, :); ccoefs(3, :); ccoefs(4, :); ];
            [i, j] = getij(3, xNe);
            obj.JT = sparse(i, j, obj.JT(:));
            % getting specefied edata
            edata = getedata(obj.etype);
            
            switch obj.etype
                case 'TTL4'
                    % evaluation of grad phi i
                    gphi = obj.JT \ repmat(edata.GG, xNe, 1);
                    % inserting in gradx and grady phi i
                    obj.gphix = (gphi(1:3:end, :))';
                    obj.gphiy = (gphi(2:3:end, :))';
                    obj.gphiz = (gphi(3:3:end, :))';
                    % evaluation of Ke
                    obj.Ke.x = zeros(10, xNe);
                    obj.Ke.y = zeros(10, xNe);
                    obj.Ke.z = zeros(10, xNe);
                    temp = 0;
                    
                    for i = 1:4
                        
                        for j = 1:i
                            temp = temp + 1;
                            obj.Ke.x(temp, :) = obj.gphix(i, :) .* obj.gphix(j, :);
                            obj.Ke.y(temp, :) = obj.gphiy(i, :) .* obj.gphiy(j, :);
                            obj.Ke.z(temp, :) = obj.gphiz(i, :) .* obj.gphiz(j, :);
                        end
                        
                    end
                    
                    obj.Ke.yx = zeros(16, xNe);
                    obj.Ke.zx = zeros(16, xNe);
                    obj.Ke.zy = zeros(16, xNe);
                    
                    for i = 1:4
                        
                        for j = 1:4
                            temp = temp + 1;
                            obj.Ke.yx(temp, :) = obj.gphiy(i, :) .* obj.gphix(j, :);
                            obj.Ke.zx(temp, :) = obj.gphiz(i, :) .* obj.gphix(j, :);
                            obj.Ke.zy(temp, :) = obj.gphiz(i, :) .* obj.gphiy(j, :);
                        end
                        
                    end
                    
                    % multiplying by triangle areas
                    obj.Ke.x = obj.Ke.x * sparse(1:xNe, 1:xNe, obj.gev);
                    obj.Ke.y = obj.Ke.y * sparse(1:xNe, 1:xNe, obj.gev);
                    obj.Ke.z = obj.Ke.z * sparse(1:xNe, 1:xNe, obj.gev);
                    obj.Ke.yx = obj.Ke.yx * sparse(1:xNe, 1:xNe, obj.gev);
                    obj.Ke.zx = obj.Ke.zx * sparse(1:xNe, 1:xNe, obj.gev);
                    obj.Ke.zy = obj.Ke.zy * sparse(1:xNe, 1:xNe, obj.gev);
                    % elemental force matrix
                    obj.Fe = repmat((obj.gev / 4), 4, 1);
                otherwise
                    error('Element type does not defined.');
            end
            
            disp('************************************')
            disp('Calculation of Ke, Fe and C completed.')
            toc
        end
        
        function addMaterial(obj, mname)
            
            if ~ ischar(mname)
                error('m name must be string.')
            end
            
            mname = lower(mname);
            mname = rmspaces(mname);
            % reading material data from binary file
            m = ReadMDataBinary(fdir, mname);
            % elaluation od axiluiry data for nonlinear magnetic materials
            if ~ m.MagneticPermeability.isLinear && m.MagneticPermeability.isIsotropic
                
                HB = m.MagneticPermeability.value;
                b = linspace(max(HB(:, 2)), 10 * max(HB(:, 2)), 100);
                b = [0.01; HB(:, 2); b(2:end)'];
                h = interp1([0; HB(:, 2)], [0; HB(:, 1)], b, 'linear', 'extrap');
                v = (h ./ b);
                
                m.BH = spline(b, h);
                m.dBdH = m.BH;
                m.dBdH.coefs = m.dBdH.coefs * diag(3:-1:1, 1);
                m.vB = spline(b, v);
                m.dvdB = m.vB;
                m.dvdB.coefs = m.dvdB.coefs * diag(3:-1:1, 1);
                
            end
            
            obj.mts.(mname) = m;
        end
        
        function setMaterial(obj, mzname, mname)
            mzname = rmspaces(mzname);
            mznames = fieldnames(obj.mzs);
            
            if ~ ismember(mzname, mznames)
                error('Specefied zone does not exist.');
            end
            
            obj.mzs.(mzname).material = lower(rmspaces(mname));
        end
        
        function set.etype(obj, etype)
            
            if ~ ischar(etype)
                error('Element type must be char.');
            end
            
            etype = upper(etype);
            
            if ~ ismember(etype, {'TTL4', 'TTL8'})
                error('Element type does not defined.');
            end
            
            obj.etype = etype;
        end
        
        %% Mesh Visiualization
        function varargout = showwf(obj)
            obj.ggmesh;
            
            f = GraphicWindow();
            f.Name = 'Wire Frame Mesh';
            h = guihandles(f);
            delete(h.bg);
            
            
            
            index = obj.facets(:, 4) ~= obj.facets(:, 5);
            patch('Faces', obj.facets(index, 1:3), 'Vertices', ...
                obj.nodes, 'FaceColor', ...
                'c', 'EdgeColor', 'w', ...
                'FaceAlpha', 0.5, 'parent', h.va);
            set(f, 'HandleVisibility', 'off', 'Visible', 'on');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
            
        end
        
        function showfb(obj, color)
            obj.ggmesh;
            
            if nargin < 2
                color = 'c';
            end
            
            ah = setFigure(['[Global Mesh Boundary Facets][', 'Nbf = ', num2str(sum(obj.bfacets)), ']'], 1);
            patch('Faces', obj.facets(obj.bfacets, 1:3), 'Vertices', obj.nodes, ...
                'FaceColor', color, 'EdgeColor', 'w', ...
                'FaceAlpha', 0.5, 'parent', ah);
            set(gcf, 'HandleVisibility', 'off', 'Visible', 'on');
        end
        
        function f = showm(obj)
            f = GraphicWindow();
            f.Name = 'Global Mesh';
            h = guihandles(f);
            mzNames = fieldnames(obj.mzs);
            
            for i = 1:numel(mzNames)
                mzptr = obj.mzs.(mzNames{i});
                mzptr.setdata;
                
                if isequal(mzptr.color, 'none')
                    edgeColor = 'none';
                else
                    edgeColor = [0.1, 0.1, 0.1];
                end
                
                patch(h.va, 'Faces', mzptr.facets(mzptr.bfacets, :), ...
                    'Vertices', mzptr.nodes, 'FaceColor', ...
                    'c', 'EdgeColor', edgeColor, ...
                    'FaceAlpha', 0.8, 'tag', mzNames{i}, 'ButtonDownFcn', @selectPatchCallback);
            end
            
            set(f, 'Visible', 'on');
        end
        
        function showmzs(obj)
            mzNames = fieldnames(obj.mzs);
            f = GraphicWindow();
            f.Name = ['[Number of Mesh Zones: ', num2str(numel(mzNames)), ']'];
            h = guihandles(f);
            
            for i = 1:numel(mzNames)
                mzptr = obj.mzs.(mzNames{i});
                mzptr.setdata;
                
                if isequal(mzptr.color, 'none')
                    edgeColor = 'none';
                else
                    edgeColor = [0.1, 0.1, 0.1];
                end
                
                patch(h.va, 'Faces', mzptr.facets(mzptr.bfacets, :), ...
                    'Vertices', mzptr.nodes, 'FaceColor', ...
                    mzptr.color, 'EdgeColor', edgeColor, ...
                    'FaceAlpha', mzptr.transparency, 'tag', mzNames{i});
            end
            
            set(f, 'Visible', 'on');
        end
        
        function showmz(obj, mzName)
            mzNames = fieldnames(obj.mzs);
            a = setFigure(['[Number of Mesh Zones: ', num2str(numel(mzNames)), ']']);
            
            for i = 1:numel(mzNames)
                mzptr = obj.mzs.(mzNames{i});
                mzptr.setdata;
                patch('Faces', mzptr.facets(mzptr.bfacets, :), ...
                    'Vertices', mzptr.nodes, 'FaceColor', ...
                    'none', 'EdgeColor', 'k', ...
                    'Parent', a, 'FaceAlpha', 0.1, 'edgeAlpha', 0.1);
            end
            
            mzptr = obj.mzs.(mzName);
            mzptr.setdata;
            patch('Faces', mzptr.facets(mzptr.bfacets, :), ...
                'Vertices', mzptr.nodes, 'FaceColor', ...
                mzptr.color, 'EdgeColor', [0.1, 0.1, 0.1], ...
                'Parent', a, 'FaceAlpha', 1);
            set(gcf, 'HandleVisibility', 'off', 'Visible', 'on');
        end
        
        function showg(obj)
            mzNames = fieldnames(obj.mzs);
            a = setFigure(['[Number of Mesh Zones: ', num2str(numel(mzNames)), ']']);
            
            for i = 1:numel(mzNames)
                mzptr = obj.mzs.(mzNames{i});
                mzptr.setdata;
                patch('Faces', mzptr.facets(mzptr.bfacets, :), ...
                    'Vertices', mzptr.nodes, 'FaceColor', ...
                    mzptr.color, 'EdgeColor', 'none', ...
                    'Parent', a, 'FaceAlpha', 1, 'FaceLighting', 'flat');
            end
            
            light('position', [1, 1, 1], 'Style', 'local')
            set(gcf, 'HandleVisibility', 'off', 'Visible', 'on');
        end
        
        function showmzss(obj)
            mznames = fieldnames(obj.mzs);
            ah = setFigure(['Number of Mesh Zones: ', num2str(numel(mznames))]);
            
            for i = 1:numel(mznames)
                if strcmpi(obj.mzs.(mznames{i}).material, 'air'), continue; end
                obj.mzs.(mznames{i}).setdata;
                patch('Faces', obj.mzs.(mznames{i}).facets(obj.mzs.(mznames{i}).bfacets, :), ...
                    'Vertices', obj.mzs.(mznames{i}).nodes, 'FaceColor', ...
                    obj.mzs.(mznames{i}).color, 'EdgeColor', 'k', ...
                    'Parent', ah, 'FaceAlpha', 1);
            end
            
            set(gca, 'HandleVisibility', 'off');
        end
        
        function showElement(obj, eindex)
            index = [obj.e(eindex, 1:4), obj.e(eindex, 6:end)];
            hold all
            
            for i = 1:length(index)
                plot3(obj.nodes(index(i), 1), obj.nodes(index(i), 2), obj.nodes(index(i), 3), ...
                    'o', 'color', 'y', 'linewidth', 1, 'MarkerFaceColor', 'y');
                text(obj.nodes(index(i), 1), obj.nodes(index(i), 2), obj.nodes(index(i), 3), num2str(i));
            end
            
            axis off equal
            
            view([1, 1, 1])
        end
        
        function varargout = shownfs(obj)
            
            f = GraphicWindow();
            f.Name = '[Named Facets]';
            h = guihandles(f);
            
            nfs = fieldnames(obj.facetNamedSelections);
            
            for i = 1:numel(nfs)
                tmp = rand(1, 3);
                patch('Faces', obj.facets(obj.facetNamedSelections.(nfs{i}), 1:3), 'Vertices', ...
                    obj.nodes, 'FaceColor', ...
                    tmp, 'EdgeColor', 'k', ...
                    'FaceAlpha', 1, 'parent', h.va);
            end
            
            legend(h.va, nfs);
            set(f, 'HandleVisibility', 'off', 'Visible', 'on');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end
        
        function varargout = shownf(obj, nfName)
            
            nfName = obj.checkFacetNamedSelectionExistence(nfName);
            
            f = GraphicWindow();
            f.Name = ['[Named Facet: ', nfName, ']'];
            h = guihandles(f);
            
            tmp = rand(1, 3);
            patch('Faces', obj.facets(obj.facetNamedSelections.(nfName), 1:3), 'Vertices', ...
                obj.nodes, 'FaceColor', ...
                tmp, 'EdgeColor', 'k', ...
                'FaceAlpha', 1, 'parent', h.va);
            legend(h.va, nfName);
            
            set(f, 'HandleVisibility', 'off', 'Visible', 'on');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end
        
        %% Named Selections
        % facet
        function name = checkFacetNamedSelectionExistence(obj, name)
            name = rmspaces(name);
            
            if ~ isfield(obj.facetNamedSelections, name)
                error('Specified facet named selection does not exist.');
            end
            
        end
        
        function name = checkFacetNamedSelectionNonExistence(obj, name)
            name = rmspaces(name);
            
            if isfield(obj.facetNamedSelections, name)
                error('Specified facet named selection already exist.');
            end
            
        end
        
        function addFacetNamedSelection(obj, name, indices)
            
            if ~ isvector(indices)
                error('indices must be a column vector.');
            end
            
            name = obj.checkFacetNamedSelectionNonExistence(name);
            obj.facetNamedSelections.(name) = indices;
            obj.facetNamedSelections.('none') = setdiff(...
                obj.facetNamedSelections.('none'), ...
                obj.facetNamedSelections.(name));
        end
        
        %% Index Finding
        function y = getfbf(obj)
            obj.ggmesh;
            y = find(obj.bfacets);
        end
        
        function y = getfbn(obj)
            y = obj.getfbf;
            y = obj.facets(y, 1:3);
            y = unique(y(:));
        end
        
        function y = getfbfiop(obj, varargin)
            y = ttmdbc_getfbfiop(obj.facets, obj.nodes, obj.facetNamedSelections.none, varargin{:});
            y = y(y ~= 0);
        end
        
        function y = getfbfiohp(obj, varargin)
            y = ttmdbc_getfbfiohp(obj.facets, obj.nodes, obj.facetNamedSelections.none, varargin{:});
            y = y(y ~= 0);
        end
        
        function y = getbfioic(obj, p0, u, r)
            y = ttmdbc_getbfioic(obj.facets, obj.nodes, obj.facetNamedSelections.none, p0, u, r);
            y = y(y ~= 0);
        end
        
        function y = getbfioc(obj, p0, u, r, zb, zu)
            y = ttmdbc_getbfioc(obj.facets, obj.nodes, obj.facetNamedSelections.none, p0, u, r, zb, zu);
            y = y(y ~= 0);
        end
        
        function y = getbfiohp(obj, p0, u1, u2)
            n = n / norm(n);
            y = ttmdbc_getbfiop(obj.facets, obj.nodes, obj.facetNamedSelections.none, p0, n);
            y = y(y ~= 0);
        end
        
        function [km, ks] = splitShift(obj, k, vec)
            Nk = length(k);
            
            if rem(Nk, 2) ~= 0
                error('number of input indices must be even.');
            end
            
            tp = obj.nodes(k, :);
            km = zeros(Nk / 2, 1);
            ks = zeros(Nk / 2, 1);
            temp = 1;
            
            while ~ isempty(tp)
                sp = tp(1, :);
                tp = tp(2:end, :);
                km(temp) = k(1);
                k = k(2:end);
                sp = sp + vec;
                index = find(sqrt(sum([tp(:, 1) - sp(1), tp(:, 2) - sp(2), tp(:, 3) - sp(3)].^2, 2)) < obj.geps);
                
                if index
                    ks(temp) = k(index);
                    tp = tp([1:index - 1, index + 1:end], :);
                    k = k([1:index - 1, index + 1:end]);
                    temp = temp + 1;
                    continue
                end
                
                sp = sp - 2 * vec;
                index = find(sqrt(sum([tp(:, 1) - sp(1), tp(:, 2) - sp(2), tp(:, 3) - sp(3)].^2, 2)) < obj.geps);
                
                if index
                    ks(temp) = km(temp);
                    km(temp) = k(index);
                    tp = tp([1:index - 1, index + 1:end], :);
                    k = k([1:index - 1, index + 1:end]);
                    temp = temp + 1;
                    continue
                end
                
                error('These set of points does not form a set of peridic points.');
            end
            
        end
        
        function y = getnIndexOnPlane(obj, p0, n)
            y = find(abs(obj.nodes * n' - p0 * n') < obj.gleps);
        end
        
        function y = getfb(obj)
            y = obj.facets(obj.bfacets, 1:3);
        end
        
        function y = geteindex(obj, k)
            
            y = ismember(obj.facets(:, 1), k) & ismember(obj.facets(:, 2), k) & ismember(obj.facets(:, 3), k);
            y = obj.facets(y, 1:3);
            
        end
        
        %% Tools
        % copy and transform
        function y = getc(obj)
            y = (obj.nodes(obj.cl(:, 1), :) + obj.nodes(obj.cl(:, 2), :) + obj.nodes(obj.cl(:, 3), :) + obj.nodes(obj.cl(:, 4), :)) / 4;
        end
        
        function copyMirrorMeshZone(obj, nmzName, mzName, varargin)
            nmzName = rmspaces(nmzName);
            mzName = rmspaces(mzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getMirror(varargin{:});
        end
        
        function copyRotateZMeshZone(obj, nmzName, mzName, varargin)
            nmzName = rmspaces(nmzName);
            mzName = rmspaces(mzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getRotateZ(varargin{:});
        end
        
        function crotatemz(obj, nmzname, mzname, varargin)
            nmzname = rmspaces(nmzname);
            mzname = rmspaces(mzname);
            obj.mzs.(nmzname) = obj.mzs.(mzname);
            obj.mzs.(nmzname).p = protatez(obj.mzs.(nmzname).p, varargin{:});
        end
        
        function cshiftmz(obj, nmzname, mzname, varargin)
            nmzname = rmspaces(nmzname);
            mzname = rmspaces(mzname);
            obj.mzs.(nmzname) = obj.mzs.(mzname).getshift(varargin{:});
        end
        
        %% Editing Functions
        function mzname = checkmzName(obj, mzname)
            mzname = rmspaces(mzname);
            
            if isfield(obj.mzs, mzname)
                error('Some mesh zones have the same name')
            end
            
        end
        
        function checkmzExistance(obj, mzname)
            mzname = rmspaces(mzname);
            
            if ~ ismember(mzname, fieldnames(obj.mzs))
                error(['Mesh zone with the name <<', mzname, ...
                    '>> does not exist.']);
            end
            
        end
        
        function y = getDefaultMeshZoneName(obj)
            index = 0;
            mzNames = fieldnames(obj.mzs);
            
            while true
                index = index + 1;
                y = ['Zone', num2str(index)];
                if ~ ismember(y, mzNames), break; end
            end
            
        end
        
        function mzname = checkMeshZoneExistence(obj, mzname)
            mzname = rmspaces(mzname);
            
            if ~ isfield(obj.mzs, mzname)
                error('Specified mesh zone does not exist.');
            end
            
        end
        
        function mzname = checkMeshZoneNonExistence(obj, mzname)
            mzname = rmspaces(mzname);
            
            if isfield(obj.mzs, mzname)
                error('Specified mesh zone already exist.');
            end
            
        end
        
        function setMeshZoneColor(obj, mzname, color)
            mzname = obj.checkMeshZoneExistence(mzname);
            obj.mzs.(rmspaces(mzname)).color = color;
        end
        
        function setmzc(varargin)
            setMeshZoneColor(varargin{:});
        end
        
        function addmz(obj, varargin)
            
            if nargin < 2, obj.EL.printError(0);
            elseif nargin == 2
                if ~ isa(varargin{1}, 'emdlab_m3d_ttmz'), obj.EL.printError(2); end
                mzname = obj.getDefaultMeshZoneName;
                mzvalue = varargin{1};
            elseif nargin == 3
                if ~ isa(varargin{1}, 'char'), obj.EL.printCharVarError('Mesh zone name'); end
                mzname = obj.checkMeshZoneNonExistence(varargin{1});
                if ~ isa(varargin{2}, 'emdlab_m3d_ttmz'), obj.EL.printError(2); end
                mzvalue = varargin{2};
            else, obj.EL.printError(1);
            end
            
            % adding new mesh zone
            obj.mzs.(mzname) = mzvalue;
            obj.mzs.(mzname).material = 'air';
            obj.mzs.(mzname).color = rand(1, 3);
            % changing states
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function removemz(obj, mzName)
            mzName = obj.checkMeshZoneExistence(mzName);
            delete(obj.mzs.(mzName));
            obj.mzs = rmfield(obj.mzs, mzName);
            % changing states
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function changeMeshZoneName(obj, mzName, newName)
            newName = obj.checkMeshZoneNonExistence(newName);
            mzName = obj.checkMeshZoneExistence(mzName);
            obj.mzs.(newName) = copy(obj.mzs.(mzName));
            obj.mzs = rmfield(obj.mzs, mzName);
        end
        
        function clearAllmzs(obj)
            obj.mzs = struct();
        end
        
        %% Topological Functions
        % setting needed data
        function obj = setdata(obj)
            % first facet of each triangle
            f1 = obj.cl(:, [1, 2, 3]);
            % second facet of each triangle
            f2 = obj.cl(:, [2, 4, 3]);
            % third facet of each triangle
            f3 = obj.cl(:, [3, 4, 1]);
            % forth facet of each triangle
            f4 = obj.cl(:, [1, 4, 2]);
            % sorting for lower index
            [f1, s1] = sort(f1, 2);
            [f2, s2] = sort(f2, 2);
            [f3, s3] = sort(f3, 2);
            [f4, s4] = sort(f4, 2);
            % specefying changed facet index
            s1 = ((s1(:, 1) == 1) & (s1(:, 2) == 3)) | ...
                ((s1(:, 1) == 3) & (s1(:, 2) == 2)) | ...
                ((s1(:, 1) == 2) & (s1(:, 2) == 1));
            s2 = ((s2(:, 1) == 1) & (s2(:, 2) == 3)) | ...
                ((s2(:, 1) == 3) & (s2(:, 2) == 2)) | ...
                ((s2(:, 1) == 2) & (s2(:, 2) == 1));
            s3 = ((s3(:, 1) == 1) & (s3(:, 2) == 3)) | ...
                ((s3(:, 1) == 3) & (s3(:, 2) == 2)) | ...
                ((s3(:, 1) == 2) & (s3(:, 2) == 1));
            s4 = ((s4(:, 1) == 1) & (s4(:, 2) == 3)) | ...
                ((s4(:, 1) == 3) & (s4(:, 2) == 2)) | ...
                ((s4(:, 1) == 2) & (s4(:, 2) == 1));
            % unification of facets
            [obj.facets, ~, ic] = unique([f1; f2; f3; f4], 'rows');
            % getting number of elements
            ne = obj.Ne;
            % getting index of facets corresponding to each elements
            f1 = ic(1:ne);
            f2 = ic(1 + ne:2 * ne);
            f3 = ic(1 + 2 * ne:3 * ne);
            f4 = ic(1 + 3 * ne:4 * ne);
            % specefying boundary facets
            obj.bfacets = sparse([f1, f2, f3, f4], ones(4 * ne, 1), ones(4 * ne, 1));
            obj.bfacets = full(obj.bfacets == 1);
            % specefying trace direction
            f1(s1) = -f1(s1);
            f2(s2) = -f2(s2);
            f3(s3) = -f3(s3);
            f4(s4) = -f4(s4);
            % element matrix
            obj.elements(:, 1:4) = [f1, f2, f3, f4];
            % edge element
            obj.facets = [obj.facets, zeros(size(obj.facets, 1), 2)];
            ttmdbc_evalfe(obj.facets, obj.elements);
            obj.facetNamedSelections.('none') = find(obj.bfacets);
        end
        
        function y = get.Nn(obj)
            y = size(obj.nodes, 1);
        end
        
        function y = get.Ne(obj)
            y = size(obj.cl, 1);
        end
        
        function y = get.Nmzs(obj)
            y = numel(fieldnames(obj.mzs));
        end
        
    end
    
    methods (Access = private)
        
        function makeFalse_isGlobalMeshGenerated(obj)
            obj.isGlobalMeshGenerated = false;
            obj.is_JIT_Evaluated = false;
            obj.isKeFe_TTL4_Evaluated = false;
            obj.isKexyz9Fe_Evaluated = false;
            obj.isKexyz1_TTL4_Evaluated = false;
            obj.isKexyz3_TTL4_Evaluated = false;
        end
        
    end
    
end
