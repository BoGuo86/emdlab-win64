% developer: https://ComProgExpert.com, Ali Jamali-Fard
% magnetic-static two-dimensional tl3
% triangular lagrangian elements: 3 points per element

classdef emdlab_solvers_ms2d_tl3 < handle
    
    properties (SetAccess = protected)
        
        % solver mesh
        m (1,1) emdlab_m2d_tmdb;

        % boundary conditions
        bcs (1,1) emdlab_bcs_scalar;

        % elements data
        edata (1,1) struct;

        % results
        results (1,1) struct;

        % excitation matrix
        exmtcs (1,1) struct;

        % dc coils
        coils (1,1) struct;
        
    end
    
    properties (SetAccess = protected)
        
        % depth of problem
        depth (1,1) double;

        % units
        units (1,1) emdlab_phy_units;

        % physical constants
        pcts (1,1) emdlab_phy_constants = emdlab_phy_constants;

        % solver Properties
        solverSettings (1,1) struct
        solverHistory (1,1) struct
        monitorResiduals (1,1) logical = true;

        % states
        isBeEvaluated (1,1) logical = false;
        isBnEvaluated (1,1) logical = false;
        isElementDataAssigned (1,1) logical = false;
        isResultsValid (1,1) logical = false;
        
    end
    
    methods        
        %% constructor and destructor
        function obj = emdlab_solvers_ms2d_tl3()
            
            % default valuess
            obj.depth = 1;
            obj.bcs = emdlab_bcs_scalar('TL3');
            obj.units = emdlab_phy_units;
            
        end
        
        function delete(obj)
            
            delete(obj.m);
            delete(obj.bcs);
            delete(obj.units);
            
        end
        
        function setUnit(obj, varargin)
            obj.units.setQuantityUnit(varargin{:});
        end
        
        function setDepth(obj, value)
            
            obj.depth = value;
            
        end
        
        %% solver properties for mesh zones
        function setdp(obj, mzName)
            
            % set default properties of a mesh zone
            obj.m.mzs.(mzName).props.isExcited = false;
            obj.m.mzs.(mzName).props.isWindingMember = false;
            obj.m.mzs.(mzName).props.isMagnetized = false;
            obj.m.mzs.(mzName).props.isCoil = false;
            obj.makeFalse_isElementDataAssigned;
            
        end
        
        function addmz(obj, mzName, mzValue)
            
            obj.m.addmz(mzName, mzValue);
            obj.setdp(mzName);
            obj.makeFalse_isElementDataAssigned;
            
        end
        
        function removemz(obj, mzName)
            
            obj.m.removemz(mzName);
            obj.makeFalse_isElementDataAssigned;
            
        end
        
        function rotateMeshZone(obj, mzName, varargin)
            
            mzName = obj.m.checkMeshZoneExistence(mzName);
            obj.m.rotateMeshZone(mzName, varargin{:});
            
            if obj.m.mzs.(mzName).props.isMagnetized
                obj.m.mzs.(mzName).props.magnetization.rotate(varargin{:});
            end
            
        end
        
        %% excitations definitions
        function windingName = checkWindingExistence(obj, windingName)
            
            windingName = rmspaces(windingName);
            
            if ~ isfield(obj.exmtcs, windingName)
                throw(MException('', 'Specified winding doen not exist.'));
            end
            
        end
        
        function windingName = checkWindingNonExistence(obj, windingName)
            
            windingName = rmspaces(windingName);
            
            if isfield(obj.exmtcs, windingName)
                throw(MException('', 'Another winding with the same name exist.'));
            end
            
        end
        
        function defineWinding(obj, windingName, varargin)
            
            windingName = obj.checkWindingNonExistence(windingName);
            obj.exmtcs.(windingName) = emdlab_solvers_ms2d_winding(varargin{:});
            % change states
            obj.makeFalse_isElementDataAssigned;
            
        end
        
        function addMeshZone2Winding(obj, windingName, mzName, mzTurns, mzDirection)
            
            windingName = obj.checkWindingExistence(windingName);
            mzName = obj.m.checkMeshZoneExistence(mzName);
            
            if obj.m.mzs.(mzName).props.isWindingMember
                throw(MException('', ['Mesh zone [', mzName, '] already is assinged to a matrix.']));
            end
            
            if obj.m.mzs.(mzName).props.isExcited
                throw(MException('', ['Mesh zone [', mzName, '] already is excited.']));
            end
            
            obj.coils.(mzName) = emdlab_solvers_ms2d_coil('turns', mzTurns, 'direction', mzDirection, ...
                'parentMatrix', windingName);
            obj.exmtcs.(windingName).addMeshZone(mzName);
            obj.m.mzs.(mzName).props.isWindingMember = true;
            % change states
            obj.makeFalse_isElementDataAssigned;
            
        end
        
        function setWindingCurrent(obj, windingName, current)
            
            windingName = obj.checkWindingExistence(windingName);
            obj.exmtcs.(windingName).current = current;
            % change states
            obj.makeFalse_isElementDataAssigned;
            
        end
        
        function setExcitation(obj, mzName, varargin)
            
            mzName = obj.m.checkMeshZoneExistence(mzName);
            
            if obj.m.mzs.(mzName).props.isWindingMember
                throw(MException('', ['Mesh zone [', mzName, '] already is assinged to a winding.']));
            end
            
            obj.m.mzs.(mzName).props.excitation = emdlab_solvers_ms2d_excitation(varargin{:});
            obj.m.mzs.(mzName).props.isExcited = true;
            % change states
            obj.makeFalse_isElementDataAssigned;
            
        end
        
        function setMagnetization(obj, mzName, varargin)
            mzName = obj.m.checkMeshZoneExistence(mzName);
            
            if nargin == 3 && isa(varargin{1}, 'emdlab_solvers_ms2d_magnetization')
                obj.m.mzs.(mzName).props.magnetization = varargin{1};
            else
                obj.m.mzs.(mzName).props.magnetization = emdlab_solvers_ms2d_magnetization(varargin{:});
            end
            
            obj.m.mzs.(mzName).props.isMagnetized = true;
            % change states
            obj.makeFalse_isElementDataAssigned;
        end
        
        %% visualization functions
        function varargout = showWinding(obj, windingName)

            f = emdlab_r2d_mesh();
            ax = axes(f);
            
            windingName = obj.checkWindingExistence(windingName);

            % winding pointer
            wptr = obj.exmtcs.(windingName);
            
            for i = 1:wptr.Nmzs
                mzptr = obj.m.mzs.(wptr.mzsName{i});
                
                if obj.coils.(wptr.mzsName{i}).sign == 1
                    mzColor = 'b';
                else
                    mzColor = 'r';
                end
                
                patch('Faces', mzptr.cl, 'Vertices', mzptr.nodes, ...
                    'FaceColor', mzColor, 'FaceAlpha', 1, 'EdgeColor', 'none', 'parent', ax);
            end
            
            index = obj.m.edges(:, 3) ~= obj.m.edges(:, 4);
            patch('Faces', obj.m.edges(index, [1, 2]), 'Vertices', obj.m.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.2, 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end
        
        function varargout = showAllWindingMembers(obj)

            f = emdlab_r2d_mesh();
            ax = axes(f);
            windingNames = fieldnames(obj.exmtcs);
            
            for i = 1:numel(windingNames)
                wptr = obj.exmtcs.(windingNames{i});
                
                for j = 1:wptr.Nmzs
                    mzptr = obj.m.mzs.(wptr.mzsName{j});
                    patch('Faces', mzptr.cl, 'Vertices', mzptr.nodes, ...
                        'FaceColor', 'g', 'FaceAlpha', 1, 'EdgeColor', 'none', 'parent', ax);
                end
                
            end

            index = obj.m.edges(:, 3) ~= obj.m.edges(:, 4);
            patch('Faces', obj.m.edges(index, [1, 2]), 'Vertices', obj.m.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.2, 'parent', ax);

            zoom on;
            axis(ax, 'off');
            axis(ax, 'equal');
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end

        end
        
        %% boundary conditions
        function clearAllBCs(obj)
            obj.bcs.clearAllBCs;
        end
        
        function setAzBC(obj, index, value, varargin)
            obj.bcs.setDirichlet(index, value, varargin{:});
        end
        
        function setOddPeriodicBC(obj, varargin)
            obj.bcs.setOddPeriodic(varargin{:});
        end
        
        function setEvenPeriodicBC(obj, varargin)
            obj.bcs.setEvenPeriodic(varargin{:});
        end
        
        %% post-proccessing: evaluations
        function clearAllResults(obj)

            % result names
            rNames = fieldnames(obj.results);
            
            for i = 1:numel(rNames)
                obj.results.(rNames{i}) = [];
            end
            
            obj.isResultsValid = false;
            
        end
        
        function evalBe(obj)

            % Evaluation of B on gaussian points
            [obj.results.Bex, obj.results.Bey] = emdlab_m2d_tl3_evalBg(obj.m.cl, obj.results.A, obj.m.JIT);
            obj.results.Bex = obj.results.Bex' * (obj.units.k_magneticVectorPotential / obj.units.k_length);
            obj.results.Bey = obj.results.Bey' * (obj.units.k_magneticVectorPotential / obj.units.k_length);

        end
        
        function evalBn(obj)
            % Evaluation of B on the mesh nodes
            if obj.isBnEvaluated, return; end
            [obj.results.Bnx, obj.results.Bny] = IHNLNRMSTL3_evalBn(obj.m.cl, obj.results.Bex, obj.results.Bey, obj.m.gea, obj.m.Nn);
            obj.isBnEvaluated = true;
        end
        
        function evalBed(obj)
            % Evaluation of B in the middle of edges
            obj.evalBn;
            obj.results.Bedx = (obj.results.Bnx(obj.m.edges(:, 1)) + obj.results.Bnx(obj.m.edges(:, 2))) / 2;
            obj.results.Bedy = (obj.results.Bny(obj.m.edges(:, 1)) + obj.results.Bny(obj.m.edges(:, 2))) / 2;
        end
                
        function [ye, yc] = evalTotalEnergyCoenergy(obj)
            
            mzsNames = fieldnames(obj.m.mzs);
            Bk2 = obj.results.Bex.^2 + obj.results.Bey.^2;
            ye = 0;
            yc = 0;
            
            for i = 1:numel(mzsNames)
                
                mzptr = obj.m.mzs.(mzsNames{i});
                mptr = obj.m.mts.(mzptr.material);
                eziptr = obj.m.ezi(:,mzptr.zi);
                
                if mptr.MagneticPermeability.isLinear
                    
                    % mesh zone energy
                    mze = 0.5*obj.edata.MagneticReluctivity(eziptr) * (mzptr.getAreaOfElements .* Bk2(eziptr));
                    ye = ye + mze;
                    yc = yc + mze;
                    
                else
                    
                    % mesh zone energy
                    mze = mzptr.getAreaOfElements' * interp1(mptr.be, mptr.we, sqrt(Bk2(eziptr)));
                    ye = ye + mze;
                    yc = yc + obj.edata.MagneticReluctivity(eziptr) * (mzptr.getAreaOfElements .* Bk2(eziptr)) - mze;
                    
                end
                
            end
            
            ye = ye*obj.getDepth*obj.units.k_length^2;
            yc = yc*obj.getDepth*obj.units.k_length^2;
            
        end
        
        function y = evalMeshZoneFluxLinkage(obj, mzName)

            mzName = obj.m.checkMeshZoneExistence(mzName);
            y = obj.m.mzs.(mzName).getQ * obj.results.A(obj.m.mzs.(mzName).l2g) * obj.getDepth;

        end
        
        function y = evalWindingFluxLinkage(obj, windingName)

            windingName = obj.checkWindingExistence(windingName);
            mptr = obj.exmtcs.(windingName);
            y = 0;
            
            for i = 1:mptr.Nmzs
                cptr = obj.coils.(mptr.mzsName{i});
                y = y + cptr.sign * cptr.turns * obj.evalMeshZoneFluxLinkage(mptr.mzsName{i});
            end
            
            y = y / mptr.np;
            
        end
        
        function y = evalTorqueByMST(obj, varargin)
            obj.evalBn;
            index = obj.m.getnIndexOnCircle(varargin{:});
            p = obj.m.nodes(index, :);
            b = [obj.results.Bnx(index), obj.results.Bny(index)];
            pAngle = atan_02pi(p) * 180 / pi;
            [~, index] = sort(pAngle);
            p = p(index, :);
            b = b(index, :);
            b = ext_xy2rt(p, b);
            b = prod(b, 2);
            b = (b + circshift(b, -1)) / 2;
            p = p - circshift(p, -1);
            p = sqrt(sum(p.^2, 2));
            y = sum(p .* b) * varargin{2} * obj.units.k_length^2 * obj.getDepth / (4 * pi * 1e-7);
        end
        
        function y = evalTorqueBySurfaceMST(obj, mzName)
            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);
            eIndex = (obj.m.edges(:, 3) == mzptr.zi) & (obj.m.edges(:, 4) ~= mzptr.zi);
            eIndex = eIndex | ((obj.m.edges(:, 4) == mzptr.zi) & (obj.m.edges(:, 3) ~= mzptr.zi));
            eIndex = eIndex & (~ obj.m.bedges);
            eIndex = find(eIndex);
            p1 = obj.m.nodes(obj.m.edges(eIndex, 1), :);
            p2 = obj.m.nodes(obj.m.edges(eIndex, 2), :);
            u = p2 - p1;
            r = (p1 + p2) / 2;
            el = sqrt(sum(u.^2, 2));
            u = [u(:, 1) ./ el, u(:, 2) ./ el];
            n = ext_protate2(u, -pi / 2);
            y = zeros(length(eIndex), 2);
            
            for i = 1:length(eIndex)
                
                if obj.m.edges(eIndex(i), 3) == mzptr.zi
                    bx = obj.results.Bex(obj.m.edges(eIndex(i), 7));
                    by = obj.results.Bey(obj.m.edges(eIndex(i), 7));
                    bt = u(i, 1) * bx + u(i, 2) * by;
                    bn = n(i, 1) * bx + n(i, 2) * by;
                    ft = bn * bt * el(i);
                    fn = 0.5 * (bn^2 - bt^2) * el(i);
                    y(i, :) = ft * u(i, :) + fn * n(i, :);
                elseif obj.m.edges(eIndex(i), 4) == mzptr.zi
                    bx = obj.results.Bex(obj.m.edges(eIndex(i), 5));
                    by = obj.results.Bey(obj.m.edges(eIndex(i), 5));
                    u(i, :) = -u(i, :);
                    n(i, :) = -n(i, :);
                    bt = u(i, 1) * bx + u(i, 2) * by;
                    bn = n(i, 1) * bx + n(i, 2) * by;
                    ft = bn * bt * el(i);
                    fn = 0.5 * (bn^2 - bt^2) * el(i);
                    y(i, :) = ft * u(i, :) + fn * n(i, :);
                else
                    error('Internal error.');
                end
                
            end
            
            y = r(:, 1) .* y(:, 2) - r(:, 2) .* y(:, 1);
            y = sum(y) * obj.units.k_length^2 * obj.getDepth / (4 * pi * 1e-7);
        end
        
        function y = evalTorqueBySurfaceMST1(obj, mzName)
            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);
            eIndex = (obj.m.edges(:, 3) == mzptr.zi) & (obj.m.edges(:, 4) ~= mzptr.zi);
            eIndex = eIndex | ((obj.m.edges(:, 4) == mzptr.zi) & (obj.m.edges(:, 3) ~= mzptr.zi));
            eIndex = eIndex & (~ obj.m.bedges);
            eIndex = find(eIndex);
            p1 = obj.m.nodes(obj.m.edges(eIndex, 1), :);
            p2 = obj.m.nodes(obj.m.edges(eIndex, 2), :);
            u = p2 - p1;
            r = (p1 + p2) / 2;
            el = sqrt(sum(u.^2, 2));
            u = [u(:, 1) ./ el, u(:, 2) ./ el];
            n = ext_protate2(u, -pi / 2);
            y = zeros(length(eIndex), 2);
            obj.evalBed;
            
            for i = 1:length(eIndex)
                
                if obj.m.edges(eIndex(i), 3) == mzptr.zi
                    bx = obj.results.Bedx(eIndex(i));
                    by = obj.results.Bedy(eIndex(i));
                    bt = u(i, 1) * bx + u(i, 2) * by;
                    bn = n(i, 1) * bx + n(i, 2) * by;
                    ft = bn * bt * el(i);
                    fn = 0.5 * (bn^2 - bt^2) * el(i);
                    y(i, :) = ft * u(i, :) + fn * n(i, :);
                elseif obj.m.edges(eIndex(i), 4) == mzptr.zi
                    bx = obj.results.Bedx(eIndex(i));
                    by = obj.results.Bedy(eIndex(i));
                    u(i, :) = -u(i, :);
                    n(i, :) = -n(i, :);
                    bt = u(i, 1) * bx + u(i, 2) * by;
                    bn = n(i, 1) * bx + n(i, 2) * by;
                    ft = bn * bt * el(i);
                    fn = 0.5 * (bn^2 - bt^2) * el(i);
                    y(i, :) = ft * u(i, :) + fn * n(i, :);
                else
                    error('Internal error.');
                end
                
            end
            
            y = r(:, 1) .* y(:, 2) - r(:, 2) .* y(:, 1);
            y = sum(y) * obj.units.k_length^2 * obj.getDepth / (4 * pi * 1e-7);
        end
        
        function y = evalTorqueByArkkio(obj, mzName, GapLength)

            obj.evalBe;
            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);
            b = [obj.results.Bex(obj.m.ezi(:, mzptr.zi)), obj.results.Bey(obj.m.ezi(:, mzptr.zi))];
            p = mzptr.getCenterOfElements;
            b = ext_xy2rt(p, b);
            r = sqrt(sum(p.^2, 2));
            y = mzptr.getAreaOfElements' * (b(:, 1) .* b(:, 2) .* r);
            y = y * obj.getDepth * obj.units.k_length^2 / GapLength / (4 * pi * 1e-7);
            
        end
        
        %% PostProccessing: Getters
        function [br, bt, tmp] = getBrBtOnCircle(obj, cx, cy, r, N)
            
            obj.evalBn;
            tmp = linspace(0, 2 * pi, N + 1);
            tmp = tmp(1:end - 1)';
            px = r * cos(tmp) + cx;
            py = r * sin(tmp) + cy;
            index = pointLocation(triangulation(obj.m.cl, obj.m.nodes), [px, py]);
            validIndex = ~ isnan(index);
            px = px(validIndex);
            py = py(validIndex);
            index = index(validIndex);
            bx = obj.results.Bex(index);
            by = obj.results.Bey(index);
            Brt = ext_xy2rt([px, py], [bx, by]);
            br = Brt(:, 1);
            bt = Brt(:, 2);
            
        end

        function [br, bt, px, py] = getBrBtOnCircleOld(obj, cx, cy, r, N)

            obj.evalBn;
            tmp = linspace(0, 2 * pi, N + 1);
            tmp = tmp(1:end - 1)';
            px = r * cos(tmp) + cx;
            py = r * sin(tmp) + cy;
            index = pointLocation(triangulation(obj.m.cl, obj.m.nodes), [px, py]);
            validIndex = ~ isnan(index);
            px = px(validIndex);
            py = py(validIndex);
            index = index(validIndex);
            Np = size(px);
            %       bx = zeros(Np, 1);
            %       by = zeros(Np, 2);
            %       for i = 1:Np
            %       end
            bx = obj.results.Bex(index);
            by = obj.results.Bey(index);
            Brt = ext_xy2rt([px, py], [bx, by]);
            br = Brt(:, 1);
            bt = Brt(:, 2);
        end
        
        %% post-proccessing: plots
        function varargout = plotAmag(obj, Ncontour)
            
            if nargin < 2
                Ncontour = 15;
            end
            
            f = figure;
            ax = axes(f);
            f.Name = 'Magnetic Vector Potential Amplitude [wb/m]';
            
            patch('faces', obj.m.cl(:, 1:3), 'vertices', obj.m.nodes, ...
                'FaceVertexCData', obj.results.A, 'FaceColor', 'interp', ...
                'EdgeColor', 'none', 'parent', ax);

            colormap(jet(Ncontour));
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Magnetic Vector Potential [A/m]';
            
            
            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch('faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', 'k', 'parent', ax);
            
            axis off equal;
            zoom on;
            set(f, 'Visible', 'on');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end

        function varargout = plotAmag3Ds(obj, Ncontour)
            
            if nargin < 2
                Ncontour = 15;
            end
            
            [f,ax] = emdlab_r3d_mesh;
            f.Name = 'Magnetic Vector Potential Amplitude [wb/m]';
            
            patch('faces', obj.m.cl(:, 1:3), 'vertices', [obj.m.nodes,obj.results.A*1000], ...
                'FaceVertexCData', obj.results.A, 'FaceColor', 'interp', ...
                'EdgeColor', 'k', 'parent', ax);

            index = obj.m.edges(:, 3) ~= obj.m.edges(:, 4);
            patch('Faces', obj.m.edges(index, [1, 2]), 'Vertices', obj.m.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.2, 'parent', ax);

            colormap(jet(Ncontour));
%             cb = colorbar(ax);
%             cb.FontName = 'Verdana';
%             cb.FontSize = 12;
%             cb.Label.String = 'Magnetic Vector Potential [A/m]';
            
            
%             index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
%             patch('faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
%                 'EdgeColor', 'k', 'parent', ax);
            
%             axis off equal;
%             zoom on;
            set(f, 'Visible', 'on');
            set(gcf,'clipping', 'off');
            set(gca,'clipping', 'off');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end
        
        function varargout = plotFluxLines(obj, N)
            
            % number of contours
            if nargin < 2
                N = 10;
            end
            
            f = figure;
            ax = axes(f);
            f.Name = 'Flux Lines';
            % evaluation of contour lines
            c = tmzpc_contour_tl3(obj.m.cl, obj.m.nodes, obj.results.A, ...
                linspace(min(obj.results.A), max(obj.results.A), N));
            t = 1:size(c, 1);
            t = reshape(t, 2, [])';
            % plotting
            index = obj.m.edges(:, 3) ~= obj.m.edges(:, 4);
            patch(ax, 'Faces', obj.m.edges(index, [1, 2]), 'Vertices', obj.m.nodes, ...
                'FaceColor', 'k', 'EdgeColor', 'k', 'LineWidth', 0.6);
            patch(ax, 'faces', t, 'vertices', c, 'edgecolor', 'b');
            set(f, 'Visible', 'on');
            axis off equal;
            zoom on;
            set(ax, 'clipping', 'off');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end
        
        function varargout = plotBmag(obj, varargin)
            
            % specefying zones
            if ~ numel(varargin)
                ti = true(obj.m.Ne, 1);
            else
                ti = false(obj.m.Ne, 1);
                
                for i = 1:numel(varargin)
                    ti = bitor(ti, obj.m.ezi(:, obj.m.mzs.(rmspaces(varargin{i})).zi));
                end
                
            end
            
            ampB = sqrt(sum(obj.results.Bex(ti).^2 + obj.results.Bey(ti).^2, 2));
            % plot through patch
            f = figure;
            ax = axes(f);
            f.Name = '[Magnetic Flux Density Amplitude [Tesla]]';
            patch(ax, 'Faces', obj.m.cl(ti, 1:3), 'Vertices', obj.m.nodes, ...
                'FaceVertexCData', ampB, 'FaceColor', 'flat', ...
                'EdgeColor', 'none');
            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', 'k');
            colormap(jet);
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Flux Density [tesla]';
%             cb.Limits(2) = 1.7;
            
            
            axis off equal;
            zoom on;
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end

        function varargout = plotBmag3D(obj, varargin)
            
            % specefying zones
            if ~ numel(varargin)
                ti = true(obj.m.Ne, 1);
            else
                ti = false(obj.m.Ne, 1);
                
                for i = 1:numel(varargin)
                    ti = bitor(ti, obj.m.ezi(:, obj.m.mzs.(rmspaces(varargin{i})).zi));
                end
                
            end
            
            ampB = sqrt(sum(obj.results.Bex(ti).^2 + obj.results.Bey(ti).^2, 2));
            % plot through patch
            f = figure;
            ax = axes(f);
            f.Name = '[Magnetic Flux Density Amplitude [Tesla]]';
            patch(ax, 'Faces', obj.m.cl(ti, 1:3), 'Vertices', [obj.m.nodes,ampB], ...
                'FaceVertexCData', ampB, 'FaceColor', 'flat', ...
                'EdgeColor', 'none');
            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', 'k');
            colormap(jet);
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Flux Density [T]';
%             cb.Limits(2) = 1.7;
            
            
            axis off equal;
            zoom on;
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end

        function varargout = plotBmagF(obj, N, varargin)
            
            % specefying zones
            if ~ numel(varargin)
                ti = true(obj.m.Ne, 1);
            else
                ti = false(obj.m.Ne, 1);
                
                for i = 1:numel(varargin)
                    ti = bitor(ti, obj.m.ezi(:, obj.m.mzs.(rmspaces(varargin{i})).zi));
                end
                
            end
            
            ampB = sqrt(sum(obj.results.Bex(ti).^2 + obj.results.Bey(ti).^2, 2));
            % plot through patch
            f = figure;
            ax = axes(f);
            f.Name = 'Magnetic Flux Density Amplitude [tesla]';
            patch(ax, 'Faces', obj.m.cl(ti, 1:3), 'Vertices', obj.m.nodes, ...
                'FaceVertexCData', ampB, 'FaceColor', 'flat', ...
                'EdgeColor', 'none');
            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', 'none');
            colormap(jet);
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Flux Density [tesla]';
            clim([0,1.9]);

            c = tmzpc_contour_tl3(obj.m.cl, obj.m.nodes, obj.results.A, ...
                linspace(min(obj.results.A), max(obj.results.A), N));
            t = 1:size(c, 1);
            t = reshape(t, 2, [])';
            patch(ax, 'faces', t, 'vertices', c, 'edgecolor', 'k');

            axis off equal;
            zoom on;
            set(ax, 'clipping', 'off');
            set(f, 'Visible', 'on');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end
        
        function varargout = plotBmagSmooth(obj, varargin)
            
            % specefying zones
            if ~ numel(varargin)
                ti = true(obj.m.Ne, 1);
            else
                ti = false(obj.m.Ne, 1);
                
                for i = 1:numel(varargin)
                    ti = bitor(ti, obj.m.ezi(:, obj.m.mzs.(rmspaces(varargin{i})).zi));
                end
                
            end
            
            f = figure;
            ax = axes(f);
            f.Name = 'Magnetic Flux Density Amplitude';
            ax.NextPlot = 'add';
            obj.evalBn;
            ampB = sqrt(obj.results.Bnx.^2 + obj.results.Bny.^2);
            % plot through patch
            patch(ax, 'Faces', obj.m.cl(ti, 1:3), 'Vertices', obj.m.nodes, ...
                'FaceVertexCData', ampB, 'FaceColor', 'interp', ...
                'EdgeColor', 'none');
            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch(ax, 'faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', 'k');
            colormap(jet);
            cb = colorbar;
            cb.FontName = 'Verdana';
            cb.FontSize = 12;
            cb.Label.String = 'Flux Density [T]';
%             cb.Limits(2) = 1.7;
            
            
            axis off equal;
            zoom on;
            set(f, 'Visible', 'on');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end
        
        function varargout = plotBvecOnCenterOfElements(obj)
            obj.evalBn;
            f = figure;
            ax = axes(f);
            f.Name = 'Magnetic field on center of mesh elements';
            ax.NextPlot = 'add';
            c = obj.m.getCenterOfElements;
            color = zeros(obj.m.Ne, 3);
            mzNames = fieldnames(obj.m.mzs);
            
            for i = 1:obj.m.Nmzs
                mzptr = obj.m.mzs.(mzNames{i});
                color(obj.m.ezi(:, mzptr.zi), :) = repmat(mzptr.color, mzptr.Ne, 1);
            end
            
            patch(ax, 'faces', obj.m.cl(:, 1:3), 'vertices', obj.m.nodes, ...
                'facecolor', 'flat', 'FaceVertexCData', color, 'EdgeColor', 'w', 'facealpha', 0.5);
            quiver(ax, c(:, 1), c(:, 2), obj.results.Bex, obj.results.Bey, 'color', 'k');
            axis(ax, 'off', 'equal');
            zoom on;
            set(f, 'Visible', 'on');
            set(ax, 'clipping', 'off');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end
        
        function varargout = plotBvecOnEdges(obj)
            
            obj.evalBed;
            f = GraphicWindow;
            f.Name = 'Magnetic field on mesh nodes';
            h = guihandles(f);
            h.va.NextPlot = 'add';
            color = zeros(obj.m.Ne, 3);
            mzNames = fieldnames(obj.m.mzs);
            
            for i = 1:obj.m.Nmzs
                mzptr = obj.m.mzs.(mzNames{i});
                color(obj.m.ezi(:, mzptr.zi), :) = repmat(mzptr.color, mzptr.Ne, 1);
            end
            
            patch(h.va, 'faces', obj.m.cl(:, 1:3), 'vertices', obj.m.nodes, ...
                'facecolor', 'flat', 'FaceVertexCData', color, 'EdgeColor', 'w', 'facealpha', 0.5);
            c = (obj.m.nodes(obj.m.edges(:, 1), :) + obj.m.nodes(obj.m.edges(:, 2), :)) / 2;
            quiver(h.va, c(:, 1), c(:, 2), obj.results.Bedx, obj.results.Bedy, 'color', 'k');
            axis(h.va, 'off', 'equal');
            set(f, 'Visible', 'on');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end
        
        function varargout = plotBvecOnNodes(obj)
            
            obj.evalBn;
            f = GraphicWindow;
            f.Name = 'Magnetic field on mesh nodes';
            h = guihandles(f);
            h.va.NextPlot = 'add';
            color = zeros(obj.m.Ne, 3);
            mzNames = fieldnames(obj.m.mzs);
            
            for i = 1:obj.m.Nmzs
                mzptr = obj.m.mzs.(mzNames{i});
                color(obj.m.ezi(:, mzptr.zi), :) = repmat(mzptr.color, mzptr.Ne, 1);
            end
            
            patch(h.va, 'faces', obj.m.cl(:, 1:3), 'vertices', obj.m.nodes, ...
                'facecolor', 'flat', 'FaceVertexCData', color, 'EdgeColor', 'w', 'facealpha', 0.5);
            quiver(h.va, obj.m.nodes(:, 1), obj.m.nodes(:, 2), obj.results.Bnx, obj.results.Bny, 'color', 'k');
            axis(h.va, 'off', 'equal');
            set(f, 'Visible', 'on');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end
        
        function plotBvecOnGrid(obj, gridSize)
            
            if nargin < 2
                gridSize = 1;
            end
            
            tmp = min(obj.m.nodes);
            xmin = tmp(1);
            ymin = tmp(2);
            tmp = max(obj.m.nodes);
            xmax = tmp(1);
            ymax = tmp(2);
            [xg, yg] = meshgrid(linspace(xmin, xmax, ceil(xmax - xmin) / gridSize), ...
                linspace(ymin, ymax, ceil(ymax - ymin) / gridSize));
            xg = xg(:);
            yg = yg(:);
            e = pointLocation(triangulation(obj.m.cl, obj.m.nodes), [xg, yg]);
            index = ~ isnan(e);
            e = e(index);
            xg = xg(index);
            yg = yg(index);
            ah = setFigure();
            quiver(xg, yg, obj.results.Bex(e), obj.results.Bey(e), 'color', 'b', 'parent', ah);
            axis off equal
            index = obj.m.edges(:, 3) - obj.m.edges(:, 4);
            patch('faces', obj.m.edges(logical(abs(index)), 1:2), 'vertices', obj.m.nodes, ...
                'EdgeColor', 'k', 'parent', ah);
            set(gcf, 'HandleVisibility', 'off', 'Visible', 'on');
        end
        
        function plotMvec(obj, scale)
            
            if nargin < 2
                scale = 1;
            end
            
            tr = triangulation(obj.m.cl(:, 1:3), obj.m.nodes);
            c = tr.incenter;
            close all;
            axis off equal; hold on;
            set(gcf, 'Color', 'k')
            quiver(c(:, 1), c(:, 2), obj.edata.MagnetizationX', obj.edata.MagnetizationY', ...
                scale, 'color', 'w');
            title('Magnetization Vector', 'color', 'c', 'fontsize', 15);
            set(gcf, 'Renderer', 'opengl');
            zoom on; hold on
        end
        
        function [p, b] = getBvecOnSegment(obj, p1, p2, N)
            if nargin < 4, N = 5; end
            % sample points
            p = zeros(N, 2);
            p(:, 1) = linspace(p1(1), p2(1), N);
            p(:, 2) = linspace(p1(2), p2(2), N);
            ti = pointLocation(triangulation(obj.m.cl, obj.m.nodes), p);
            b = [obj.results.Bex(ti), obj.results.Bey(ti)];
        end
        
        function [xp, xB] = getBOnCircle(obj, c, r, N)
            tr = triangulation(obj.m.cl(:, 1:3), obj.m.nodes);
            xp = zeros(N, 2);
            
            for i = 1:N
                xp(i, :) = protate(c + [r, 0], (i - 1) * 2 * pi / N, c);
            end
            
            index = tr.pointLocation(xp);
            index = index(~ isnan(index));
            xB = obj.B(index, :);
        end
        
        function plotBvecOnCircle(obj, c, r, N)
            [xp, xB] = obj.getBrBtOnCircleOld(c, r, N);
            close all
            quiver(xp(:, 1), xp(:, 2), xB(:, 1), xB(:, 2));
            axis off equal
        end
        
        function plotBrBtOnCircle(obj, varargin)
            %       obj.evalBn;
            %       index = obj.m.getnIndexOnCircle(varargin{:});
            %       p = obj.m.nodes(index,:);
            %       b = [obj.results.Bnx(index),obj.results.Bny(index)];
            %       pAngle = atan_02pi(p)*180/pi;
            %       [pAngle,index] = sort(pAngle);
            %       p = p(index,:);
            %       b = b(index,:);
            %       b = ext_xy2rt(p,b);
            [br, bt, ~, ~] = obj.getBrBtOnCircleOld(varargin{:});
            figure('Name', 'Br and Bt on circle')
            subplot(211)
            t = linspace(0,360,length(br));
            plot(t,br)
                  set(gca,'Xlim',[0,360])
            xlabel('Mechanical Angle [Degree]')
            ylabel('B_r [Tesla]')
            title('Rdaial Flux Density Waveform');
            subplot(212)
            plot(t,bt)
                  set(gca,'Xlim',[0,360])
            xlabel('Mechanical Angle [Degree]')
            ylabel('B_t [Tesla]')
            title('Tangentional Flux Density Waveform');
        end
        
        function plotBrfft(obj, POLES, Nh, varargin)
            obj.evalBn;
            index = obj.m.getnIndexOnCircle(varargin{:});
            p = obj.m.nodes(index, :);
            b = [obj.Bn.x(index), obj.Bn.y(index)];
            pAngle = atan_02pi(p);
            [pAngle, index] = sort(pAngle);
            p = p(index, :);
            b = b(index, :);
            b = ext_xy2rt(p, b);
            figure('Name', '[@EMDLab] Middle Air Gap Br and Bt', 'NumberTitle', 'off', ...
                'WindowStyle', 'modal')
            subplot(211)
            plot(pAngle * 180 / pi, b(:, 1))
            set(gca, 'Xlim', [pAngle(1) * 180 / pi, pAngle(end) * 180 / pi])
            xlabel('Mechanical Angle [Degree]')
            ylabel('B_r [Tesla]')
            title(['Rdaial Flux Density, Mean = ', num2str(mean(abs(b(:, 1))))])
            hold all;
            an = (1 / pi) * int_trap(full(repmat(b(:, 1)', Nh, 1) .* cos(POLES * sparse(1:Nh, 1:Nh, 1:Nh) * repmat(pAngle', Nh, 1))), pAngle');
            bn = (1 / pi) * int_trap(full(repmat(b(:, 1)', Nh, 1) .* sin(POLES * sparse(1:Nh, 1:Nh, 1:Nh) * repmat(pAngle', Nh, 1))), pAngle');
            hn = sqrt(an.^2 + bn.^2);
            y = zeros(length(pAngle), 1);
            
            for i = 1:Nh
                y = y + an(i) * cos(POLES * i * pAngle) + bn(i) * sin(POLES * i * pAngle);
                plot(pAngle * 180 / pi, y, '-.');
            end
            
            subplot(212)
            stem(1:Nh, hn)
            set(gca, 'Xlim', [0, Nh + 1])
            xlabel('Harmonic Index')
            ylabel('Harmonic Amplitude')
            title(['Harmonic Decomposition of Radial Flux, H[1] = ', num2str(hn(1))]);
        end
        
        function y = getbonb(obj, mzname, bname)
            b = obj.m.mzs.(mzname).bs.(bname);
            
            p = [b.p1; b.getip; b.p2];
            
            tr = triangulation(obj.m.cl(:, 1:3), obj.m.nodes);
            [~, index] = ismembertol(p, obj.m.nodes, 1e-6, 'ByRows', true);
            y = tr.vertexAttachments(index);
            b = zeros(numel(y), 2);
            
            for i = 1:numel(y)
                tmp = 0;
                
                for j = 1:length(y{i})
                    tmp = tmp + obj.B(y{i}(j), :) * obj.m.gea(y{i}(j));
                end
                
                b(i, :) = tmp / sum(obj.m.gea(y{i}));
            end
            
            y = b;
        end
        
        %% PostProccessing: Error estimator
        function y = getElementsError(obj)
            y = 0;
        end
        
    end
    
    % solver methods
    methods
        
        function y = getTotalEnergy(obj)
            y = obj.evalTotalEnergy;
        end
        
        function y = getDepth(obj)
            y = obj.depth * obj.units.k_length;
        end
        
        function y = getMeshZoneAc(obj, mzName)
            mzName = obj.m.checkMeshZoneExistence(mzName);
            mzptr = obj.m.mzs.(mzName);
            y = obj.results.A(mzptr.l2g);
            y = y(mzptr.cl);
            y = mean(y, 2);
        end
        
    end
    
    methods (Access = protected)
        
        function makeFalse_isElementDataAssigned(obj)
            obj.isElementDataAssigned = false;
            obj.isBnEvaluated = false;
        end
        
    end
    
end
