% developer: https://ComProgExpert.com, Ali Jamali-Fard
% 2D triangular mesh data base

classdef emdlab_m2d_tmdb < handle & emdlab_g2d_constants & matlab.mixin.Copyable & emdlab_mdb_cp
    
    properties (SetAccess = private)
        
        % mesh nodes
        nodes (:,2) double;

        % mesh connectivity list
        cl (:,:) double;

        % mesh elements: [edge1, edge2, edge3, zone index]
        elements (:,4) double;

        % unique edges (:,8): [index of node1, index of node2, ]
        edges

        % list of boundary edges
        bedges

        % jacobian inverse transpose
        JIT (4,:) double;

        % global element area
        gea (1,:) double;

        % element zone index
        ezi (:,:) logical;
        
        % auxiliary stored matricies
        mtcs (1,1) struct;

        % named selections
        edgeNamedSelections (1,1) struct;

        % flag to print the elapsed times
        printFlag (1,1) logical = true;
        
    end
    
    properties (Dependent = true)
        
        % Number of nodes
        Nn (1,1) double;
        
        % Number of elements
        Ne (1,1) double;
        
    end
    
    properties (Access = private)
        
        % element type
        etype (1,:) char;
        isd2ElementsGenerated (1,1) logical = false;
        isd3ElementsGenerated (1,1) logical = false;
        isJITEvaluated (1,1) logical = false;
        isKeMeFe_TL3_Evaluated (1,1) logical = false;
        isKeMeFe_TL6_Evaluated (1,1) logical = false;
        isKeFe_TL3_Evaluated (1,1) logical = false;
        isKexy4Fe_TL3_Evaluated (1,1) logical = false;
        isKeFe_TL6_Evaluated (1,1) logical = false;
        isKeFe_TL6_cte_Evaluated (1,1) logical = false;
        isKe_TL3_Evaluated (1,1) logical = false;
        isFe_TL3_Evaluated (1,1) logical = false;
        isMe_TL3_Evaluated (1,1) logical = false;
        isKexy1_TL3_Evaluated (1,1) logical = false;
        isKexy4_TL3_Evaluated (1,1) logical = false;
        isKexy1_TL6_Evaluated (1,1) logical = false;
        isKerz1_TL3_Evaluated (1,1) logical = false;
        
    end
    
    methods
        %% Constructor and Destructor
        function obj = emdlab_m2d_tmdb()
        end
        
        function y = get.Nn(obj)
            y = size(obj.nodes, 1);
        end
        
        function y = get.Ne(obj)
            y = size(obj.cl, 1);
        end
        
        function setPrintFlag(obj, newValue)
            obj.printFlag = newValue;
        end

        function delete(obj)
            
            mzNames = fieldnames(obj.mzs);
            
            for i = 1:numel(mzNames)
                delete(obj.mzs.(mzNames{i}));
            end
            
        end
        
        %% FEM Preparation
        function ggmesh(obj)
            % generate|global|mesh
            
            % check states
            if obj.isGlobalMeshGenerated, return; end
            
            % generation of initial mesh
            Nn_tmp = 0;
            Ne_tmp = 0;
            mzNames = fieldnames(obj.mzs);
            
            % loop over mesh zones: for calculation number of total nodes and elements
            for i = 1:obj.Nmzs
                Nn_tmp = Nn_tmp + obj.mzs.(mzNames{i}).Nn;
                Ne_tmp = Ne_tmp + obj.mzs.(mzNames{i}).Ne;
            end
            
            % initialization of nodes and elements
            obj.nodes = zeros(Nn_tmp, 2);
            obj.cl = zeros(Ne_tmp, 3);
            obj.elements = zeros(Ne_tmp, 4);
            nindex = 0;
            eindex = 0;
            
            % loop over mesh zones: for insertion of nodes and elements
            for i = 1:obj.Nmzs
                
                % insertion of nodes
                obj.nodes(1 + nindex:obj.mzs.(mzNames{i}).Nn + nindex, :) = ...
                    obj.mzs.(mzNames{i}).nodes;
                
                % insertion of elements
                obj.cl(1 + eindex:obj.mzs.(mzNames{i}).Ne + eindex, 1:3) = ...
                    obj.mzs.(mzNames{i}).cl + nindex;
                
                % specefying zone index
                obj.elements(1 + eindex:obj.mzs.(mzNames{i}).Ne + eindex, 4) = i;
                obj.mzs.(mzNames{i}).zi = i;
                nindex = nindex + obj.mzs.(mzNames{i}).Nn;
                eindex = eindex + obj.mzs.(mzNames{i}).Ne;
                
            end
            
            % unification of nodes
            [obj.nodes, ~, ic] = uniquetol(obj.nodes, obj.gleps, 'ByRows', true);
            obj.cl = ic(obj.cl);
            
            % setting l2g
            nindex = 0;
            
            % loop over mesh zones: for setting l2g
            for i = 1:obj.Nmzs
                obj.mzs.(mzNames{i}).l2g = ic(nindex + 1:nindex + ...
                    obj.mzs.(mzNames{i}).Nn);
                nindex = nindex + obj.mzs.(mzNames{i}).Nn;
            end
            
            obj.setdata;
            obj.evalezi;
            
            % change states
            obj.isGlobalMeshGenerated = true;
            
            % settig element type
            obj.etype = 'TL3';
            
        end
        
        % evaluate element zone index
        function evalezi(obj)
            
            obj.ezi = false(obj.Ne, obj.Nmzs);
            
            for i = 1:obj.Nmzs
                obj.ezi(:, i) = obj.elements(:, 4) == i;
            end
            
        end
        
        % evaluate jacobian inverse transpose matrix
        function evalJIT(obj)
                        
            % check states
            if obj.isJITEvaluated, return; end
            if obj.printFlag
                tic, disp('-------------------------------------------------------');
            end
            
            % prerequisite
            obj.ggmesh;
            
            % evaluation of jacobian inverse transpose using d1 element data
            % JIT is a [4 x Ne] matrix [J11;J21;J12;J22]
            [obj.JIT, obj.gea] = emdlab_m2d_tl3_evalJIT(obj.cl(:,1:3), obj.nodes);
            
            % change states
            obj.isJITEvaluated = true;
            if obj.printFlag
                disp('Evaluation of JIT completed.');
                toc, disp('-------------------------------------------------------');
            end
            
        end

        function evalJIT_old(obj)
                        
            % check states
            if obj.isJITEvaluated, return; end
            if obj.printFlag
                tic, disp('-------------------------------------------------------');
            end
            
            % prerequisite
            obj.ggmesh;
            
            % gea: global element vector
            obj.gea = zeros(1, obj.Ne);
            mzNames = fieldnames(obj.mzs);
            
            % loop over mesh zones: for gea: global element area
            for i = 1:obj.Nmzs
                obj.gea(obj.ezi(:, obj.mzs.(mzNames{i}).zi)) = ...
                    obj.mzs.(mzNames{i}).getAreaOfElements;
            end
            
            % evaluation of jacobian inverse transpose using d1 element data
            % JIT is a [4xNe] matrix [J11;J21;J12;J22]
            obj.JIT = tmdbc_evalJIT(obj.cl, obj.nodes, obj.gea);
            
            % change states
            obj.isJITEvaluated = true;
            if obj.printFlag
                disp('Evaluation of JIT completed.');
                toc, disp('-------------------------------------------------------');
            end
            
        end
        
        % triangular elements: first order
        function evalKeMeFe_TL3(obj)

            % prerequests
            obj.ggmesh;
            obj.evalJIT;
            
            tic, disp('-------------------------------------------------------');
            if obj.isKeFe_TL3_Evaluated, return; end

            % calculation of stiffness matrix, mass matrix and force vector
            [obj.mtcs.Ke, obj.mtcs.Me, obj.mtcs.Fe] = emdlab_m2d_tl3_evalKeMeFe(obj.cl, obj.nodes);
            
            disp('Calculation of [Ke], [Me], and [Fe] completed.');
            toc, disp('-------------------------------------------------------');
            
            % change states
            obj.isKeMeFe_TL3_Evaluated = true;

        end

        function evalKeMeFe_TL6(obj)

            % prerequests
            obj.ggmesh;
            obj.evalJIT;
            
            tic, disp('-------------------------------------------------------');
            if obj.isKeFe_TL6_Evaluated, return; end

            % calculation of stiffness matrix, mass matrix and force vector
            [obj.mtcs.Ke1, obj.mtcs.Ke2, obj.mtcs.Ke3, obj.mtcs.Me, obj.mtcs.Fe, obj.mtcs.FeMx, obj.mtcs.FeMy] = emdlab_m2d_tl6_evalKeMeFe(obj.cl, obj.nodes);
            
            disp('Calculation of [Ke], [Me], and [Fe] completed.');
            toc, disp('-------------------------------------------------------');
            
            % change states
            obj.isKeMeFe_TL6_Evaluated = true;

        end

        function evalKeFe_TL3(obj)
            
            % prerequests
            obj.ggmesh;
            obj.evalJIT;
            
            tic, disp('-------------------------------------------------------');
            if obj.isKeFe_TL3_Evaluated, return; end
            
            % getting number of elements
            xNe = obj.Ne;
            
            % evaluation of grad phi i
            [obj.gphix, obj.gphiy] = tmdbc_evalGphixy_TL3(obj.JIT);
            
            % evaluation of Ke
            obj.Ke = zeros(6, xNe);
            temp = 0;
            
            for i = 1:3
                
                for j = 1:i
                    temp = temp + 1;
                    obj.Ke(temp, :) = obj.gphix(i, :) .* obj.gphix(j, :) + ...
                        obj.gphiy(i, :) .* obj.gphiy(j, :);
                end
                
            end
            
            % multiplying by triangle areas
            obj.Ke = obj.Ke * sparse(1:xNe, 1:xNe, obj.gea);
            
            % elemental force matrix
            obj.Fe = repmat((obj.gea / 3), 3, 1);
            disp('Calculation of [Ke] and [Fe] completed.');
            toc, disp('-------------------------------------------------------');
            
            % change states
            obj.isKeFe_TL3_Evaluated = true;

        end
        
        function evalKexy1_TL3(obj)
            
            if obj.isKexy1_TL3_Evaluated, return; end
            
            obj.ggmesh;
            obj.evalJIT;
            if obj.printFlag
                tic, disp('-------------------------------------------------------');
            end
            
            % evaluation of grad phi i
            [obj.mtcs.gphix, obj.mtcs.gphiy] = tmdbc_evalGphixy_TL3(obj.JIT);
            
            % evaluation of Ke
            obj.mtcs.Ke = zeros(6, obj.Ne);
            tmp = 0;
            
            for i = 1:3
                
                for j = 1:i
                    tmp = tmp + 1;
                    obj.mtcs.Ke(tmp, :) = (obj.mtcs.gphix(i, :) .* obj.mtcs.gphix(j, :) + ...
                        obj.mtcs.gphiy(i, :) .* obj.mtcs.gphiy(j, :)) .* obj.gea;
                end
                
            end
            
            if obj.printFlag
                disp('Calculation of [Kexy1] completed.');
                toc, disp('-------------------------------------------------------');
            end

            % change states
            obj.isKexy1_TL3_Evaluated = true;
            
        end
        
        function evalKerz1_TL3(obj)
            if obj.isKerz1_TL3_Evaluated, return; end
            obj.ggmesh;
            obj.evalJIT;
            tic, disp('-------------------------------------------------------');
            % evaluation of grad phi i
            [obj.mtcs.gphix, obj.mtcs.gphiy] = tmdbc_evalGphixy_TL3(obj.JIT);
            % evaluation of Ke
            obj.mtcs.Ke = zeros(6, obj.Ne);
            tmp = 0;
            
            for i = 1:3
                
                for j = 1:i
                    tmp = tmp + 1;
                    obj.mtcs.Ke(tmp, :) = (obj.mtcs.gphix(i, :) .* obj.mtcs.gphix(j, :) + ...
                        obj.mtcs.gphiy(i, :) .* obj.mtcs.gphiy(j, :)) .* obj.gea;
                end
                
            end
            
            disp('Calculation of [Kexy1] completed.');
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKerz1_TL3_Evaluated = true;
        end
        
        function evalKexy4_TL3(obj)
            if obj.isKexy4_TL3_Evaluated, return; end
            obj.ggmesh;
            obj.evalJIT;
            tic, disp('-------------------------------------------------------');
            % evaluation of grad phi i
            [obj.mtcs.gphix, obj.mtcs.gphiy] = tmdbc_evalGphixy_TL3(obj.JIT);
            % evaluation of Ke
            obj.mtcs.Kex = zeros(6, obj.Ne);
            obj.mtcs.Key = zeros(6, obj.Ne);
            obj.mtcs.Kexy = zeros(9, obj.Ne);
            obj.mtcs.Keyx = zeros(9, obj.Ne);
            tmp = 0;
            
            for i = 1:3
                
                for j = 1:i
                    tmp = tmp + 1;
                    obj.mtcs.Kex(tmp, :) = obj.mtcs.gphix(i, :) .* obj.mtcs.gphix(j, :) .* obj.gea;
                    obj.mtcs.Key(tmp, :) = obj.mtcs.gphiy(i, :) .* obj.mtcs.gphiy(j, :) .* obj.gea;
                end
                
            end
            
            tmp = 0;
            
            for i = 1:3
                
                for j = 1:3
                    tmp = tmp + 1;
                    obj.mtcs.Kexy(tmp, :) = obj.mtcs.gphix(i, :) .* obj.mtcs.gphiy(j, :) .* obj.gea;
                    obj.mtcs.Keyx(tmp, :) = obj.mtcs.gphiy(i, :) .* obj.mtcs.gphix(j, :) .* obj.gea;
                end
                
            end
            
            disp('Calculation of [Kexy4] completed.');
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKexy4_TL3_Evaluated = true;
        end
        
        function evalKexy1_TL6_isoparametric(obj)
            
            if obj.isKexy1_TL6_Evaluated, return; end
            
            obj.ggmesh;
            obj.gd2elements;
            obj.evalJIT;
            tic, disp('-------------------------------------------------------');
            
            % getting number of elements
            xNe = obj.Ne;
            
            % getting data of reference element
            tic
            % getting specefied edata
            edata = getedata('TL6');
            
            % JITGG at six point of TR
            % point1
            tmp = obj.JIT * repmat(edata.GG(0, 0), xNe, 1);
            obj.mtcs.JITGG00x = tmp(1:2:end, :)';
            obj.mtcs.JITGG00y = tmp(2:2:end, :)';
            % point 2
            tmp = obj.JIT * repmat(edata.GG(1, 0), xNe, 1);
            obj.mtcs.JITGG10x = tmp(1:2:end, :)';
            obj.mtcs.JITGG10y = tmp(2:2:end, :)';
            % point 3
            tmp = obj.JIT * repmat(edata.GG(0, 1), xNe, 1);
            obj.mtcs.JITGG01x = tmp(1:2:end, :)';
            obj.mtcs.JITGG01y = tmp(2:2:end, :)';
            % point 4
            tmp = obj.JIT * repmat(edata.GG(0.5, 0), xNe, 1);
            obj.mtcs.JITGG120x = tmp(1:2:end, :)';
            obj.mtcs.JITGG120y = tmp(2:2:end, :)';
            % point 5
            tmp = obj.JIT * repmat(edata.GG(0.5, 0.5), xNe, 1);
            obj.mtcs.JITGG1212x = tmp(1:2:end, :)';
            obj.mtcs.JITGG1212y = tmp(2:2:end, :)';
            % point 6
            tmp = obj.JIT * repmat(edata.GG(0, 0.5), xNe, 1);
            obj.mtcs.JITGG012x = tmp(1:2:end, :)';
            obj.mtcs.JITGG012y = tmp(2:2:end, :)';
            % evaluation of Ke
            obj.Ke = zeros(21, xNe);
            tmp = 0;
            
            for i = 1:6
                
                for j = 1:i
                    tmp = tmp + 1;
                    obj.Ke(tmp, :) = (obj.mtcs.JITGG120x(j, :)) .* (obj.mtcs.JITGG120x(i, :)) + ...
                        (obj.mtcs.JITGG120y(j, :)) .* (obj.mtcs.JITGG120y(i, :)) + ...
                        (obj.mtcs.JITGG1212x(j, :)) .* (obj.mtcs.JITGG1212x(i, :)) + ...
                        (obj.mtcs.JITGG1212y(j, :)) .* (obj.mtcs.JITGG1212y(i, :)) + ...
                        (obj.mtcs.JITGG012x(j, :)) .* (obj.mtcs.JITGG012x(i, :)) + ...
                        (obj.mtcs.JITGG012y(j, :)) .* (obj.mtcs.JITGG012y(i, :));
                end
                
            end
            
            % multiplying by triangle areas
            obj.Ke = obj.Ke * sparse(1:xNe, 1:xNe, obj.gea / 3);
            % calculation of fi
            obj.Fe = [zeros(3, xNe); repmat((obj.gea / 3), 3, 1)];
            disp('Calculation of [Ke] completed.')
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKexy1_TL6_Evaluated = true;
        end
        
        function evalKexy1_TL6(obj)
            
            % prerequests
            obj.ggmesh;
            obj.evalJIT;
            
            tic, disp('-------------------------------------------------------');
            
            if obj.isKeFe_TL6_cte_Evaluated, return; end
            
            % getting number of elements
            xNe = obj.Ne;
            % getting data of reference element
            tic
            % getting specefied edata
            edata = getedata('TL6');
            % GG at three point of TR
            % at first point
            obj.mtcs.GG120 = repmat(edata.GG(0.5, 0), xNe, 1);
            % at second point
            obj.mtcs.GG1212 = repmat(edata.GG(0.5, 0.5), xNe, 1);
            % at third point
            obj.mtcs.GG012 = repmat(edata.GG(0, 0.5), xNe, 1);
            % evaluation of Ke
            obj.mtcs.Ke = zeros(21, xNe);
            temp = 0;
            
            for i = 1:6
                
                for j = i:6
                    temp = temp + 1;
                    ke = (obj.JIT * obj.mtcs.GG120(:, j)) .* (obj.JIT * obj.mtcs.GG120(:, i)) + ...
                        (obj.JIT * obj.mtcs.GG1212(:, j)) .* (obj.JIT * obj.mtcs.GG1212(:, i)) + ...
                        (obj.JIT * obj.mtcs.GG012(:, j)) .* (obj.JIT * obj.mtcs.GG012(:, i));
                    ke = reshape(ke, 2, []);
                    obj.mtcs.Ke(temp, :) = sum(ke);
                end
                
            end
            
            % multiplying by triangle areas
            obj.mtcs.Ke = obj.mtcs.Ke * sparse(1:xNe, 1:xNe, obj.gea / 3);
            % calculation of fi
            obj.mtcs.Fe = [zeros(3, xNe); repmat((obj.gea / 3), 3, 1)];
            disp('Calculation of [Ke] and [Fe] completed.');
            toc, disp('-------------------------------------------------------');
            % change states
            obj.isKexy1_TL6_Evaluated = true;
        end
        
        function evalKeFe_TL6_cte(obj)
            obj.ggmesh;
            evalJIT;
            
            if obj.isKeFe_TL6_cte_Evaluated
                return
            end
            
            % getting number of elements
            xNe = obj.Ne;
            % getting data of reference element
            tic
            % getting specefied edata
            edata = getedata('TL6');
            % GG at three point of TR
            % at first point
            GG120 = repmat(edata.GG(0.5, 0), xNe, 1);
            % at second point
            GG1212 = repmat(edata.GG(0.5, 0.5), xNe, 1);
            % at third point
            GG012 = repmat(edata.GG(0, 0.5), xNe, 1);
            % evaluation of Ke
            obj.Ke = zeros(21, xNe);
            temp = 0;
            
            for i = 1:6
                
                for j = i:6
                    temp = temp + 1;
                    ke = (obj.JIT * GG120(:, j)) .* (obj.JIT * GG120(:, i)) + ...
                        (obj.JIT * GG1212(:, j)) .* (obj.JIT * GG1212(:, i)) + ...
                        (obj.JIT * GG012(:, j)) .* (obj.JIT * GG012(:, i));
                    ke = reshape(ke, 2, []);
                    obj.Ke(temp, :) = sum(ke);
                end
                
            end
            
            % multiplying by triangle areas
            obj.Ke = obj.Ke * sparse(1:xNe, 1:xNe, obj.gea / 3);
            % calculation of fi
            obj.Fe = [zeros(3, xNe); repmat((obj.gea / 3), 3, 1)];
            disp('************************************')
            disp('Calculation of Ke, Fe and C completed.')
            toc
            % change states
            obj.isKeFe_TL6_cte_Evaluated = true;
        end
        
        function evalKeFe_TL6(obj)
            if obj.isKeFe_TL6_Evaluated, return; end
            obj.ggmesh;
            obj.evalJIT;
            % getting number of elements
            xNe = obj.Ne;
            % getting data of reference element
            tic
            % getting specefied edata
            edata = getedata('TL6');
            % JITGG at six point of TR
            % point1
            tmp = obj.JIT * repmat(edata.GG(0, 0), xNe, 1);
            obj.mtcs.JITGG00x = tmp(1:2:end, :)';
            obj.mtcs.JITGG00y = tmp(2:2:end, :)';
            % point 2
            tmp = obj.JIT * repmat(edata.GG(1, 0), xNe, 1);
            obj.mtcs.JITGG10x = tmp(1:2:end, :)';
            obj.mtcs.JITGG10y = tmp(2:2:end, :)';
            % point 3
            tmp = obj.JIT * repmat(edata.GG(0, 1), xNe, 1);
            obj.mtcs.JITGG01x = tmp(1:2:end, :)';
            obj.mtcs.JITGG01y = tmp(2:2:end, :)';
            % point 4
            tmp = obj.JIT * repmat(edata.GG(0.5, 0), xNe, 1);
            obj.mtcs.JITGG120x = tmp(1:2:end, :)';
            obj.mtcs.JITGG120y = tmp(2:2:end, :)';
            % point 5
            tmp = obj.JIT * repmat(edata.GG(0.5, 0.5), xNe, 1);
            obj.mtcs.JITGG1212x = tmp(1:2:end, :)';
            obj.mtcs.JITGG1212y = tmp(2:2:end, :)';
            % point 6
            tmp = obj.JIT * repmat(edata.GG(0, 0.5), xNe, 1);
            obj.mtcs.JITGG012x = tmp(1:2:end, :)';
            obj.mtcs.JITGG012y = tmp(2:2:end, :)';
            % evaluation of Ke
            obj.Ke = zeros(21, xNe);
            tmp = 0;
            
            for i = 1:6
                
                for j = 1:i
                    tmp = tmp + 1;
                    obj.Ke(tmp, :) = (obj.mtcs.JITGG120x(j, :)) .* (obj.mtcs.JITGG120x(i, :)) + ...
                        (obj.mtcs.JITGG120y(j, :)) .* (obj.mtcs.JITGG120y(i, :)) + ...
                        (obj.mtcs.JITGG1212x(j, :)) .* (obj.mtcs.JITGG1212x(i, :)) + ...
                        (obj.mtcs.JITGG1212y(j, :)) .* (obj.mtcs.JITGG1212y(i, :)) + ...
                        (obj.mtcs.JITGG012x(j, :)) .* (obj.mtcs.JITGG012x(i, :)) + ...
                        (obj.mtcs.JITGG012y(j, :)) .* (obj.mtcs.JITGG012y(i, :));
                end
                
            end
            
            % multiplying by triangle areas
            obj.Ke = obj.Ke * sparse(1:xNe, 1:xNe, obj.gea / 3);
            % calculation of fi
            obj.Fe = [zeros(3, xNe); repmat((obj.gea / 3), 3, 1)];
            disp('************************************')
            disp('Calculation of Ke, Fe and C completed.')
            toc
            % change states
            obj.isKeFe_TL6_Evaluated = true;
        end
        
        function evalKeFe(obj, etype)
            obj.etype = etype;
            % getting number of elements
            xNe = obj.Ne;
            % getting data of reference element
            tic
            % evaluation of jacobian inverse using d1 element data
            edata = getedata('TL3');
            % x and y coordinate of points
            xp = obj.nodes(:, 1);
            yp = obj.nodes(:, 2);
            % point coordinate of each triangle nodes
            xp = xp(obj.cl(:, 1:3)');
            yp = yp(obj.cl(:, 1:3)');
            % calculation of area of each triangles
            obj.gea = zeros(1, xNe);
            mznames = fieldnames(obj.mzs);
            
            for i = 1:obj.Nmzs
                obj.gea(obj.ezi(:, obj.mzs.(mznames{i}).zi)) = ...
                    obj.mzs.(mznames{i}).getAreaOfElements;
            end
            
            % evaluation of a and b
            acoefs = edata.M \ xp;
            bcoefs = edata.M \ yp;
            % Jacobian inverse transpose of each elements
            obj.JIT = [bcoefs(3, :); -acoefs(3, :); ...
                -bcoefs(2, :); acoefs(2, :)];
            obj.JIT = obj.JIT * sparse(1:xNe, 1:xNe, 1 ./ (2 * obj.gea));
            [i, j] = getij(2, xNe);
            obj.JIT = sparse(i, j, obj.JIT(:));
            % getting specefied edata
            edata = getedata(obj.etype);
            
            switch upper(obj.etype)
                case 'TL3'
                    % evaluation of grad phi i
                    gphi = obj.JIT * repmat(edata.GG, xNe, 1);
                    % inserting in gradx and grady phi i
                    obj.gphix = (gphi(1:2:end, :))';
                    obj.gphiy = (gphi(2:2:end, :))';
                    % evaluation of Ke
                    obj.Ke = zeros(6, xNe);
                    temp = 0;
                    
                    for i = 1:3
                        
                        for j = 1:i
                            temp = temp + 1;
                            obj.Ke(temp, :) = obj.gphix(i, :) .* obj.gphix(j, :) + ...
                                obj.gphiy(i, :) .* obj.gphiy(j, :);
                        end
                        
                    end
                    
                    % multiplying by triangle areas
                    obj.Ke = obj.Ke * sparse(1:xNe, 1:xNe, obj.gea);
                    % elemental force matrix
                    obj.Fe = repmat((obj.gea / 3), 3, 1);
                case 'TL6'
                    % GG at three point of TR
                    % at first point
                    GG120 = repmat(edata.GG(0.5, 0), xNe, 1);
                    % at second point
                    GG1212 = repmat(edata.GG(0.5, 0.5), xNe, 1);
                    % at third point
                    GG012 = repmat(edata.GG(0, 0.5), xNe, 1);
                    % evaluation of Ke
                    obj.Ke = zeros(21, xNe);
                    temp = 0;
                    
                    for i = 1:6
                        
                        for j = i:6
                            temp = temp + 1;
                            ke = (obj.JIT * GG120(:, j)) .* (obj.JIT * GG120(:, i)) + ...
                                (obj.JIT * GG1212(:, j)) .* (obj.JIT * GG1212(:, i)) + ...
                                (obj.JIT * GG012(:, j)) .* (obj.JIT * GG012(:, i));
                            ke = reshape(ke, 2, []);
                            obj.Ke(temp, :) = sum(ke);
                        end
                        
                    end
                    
                    % multiplying by triangle areas
                    obj.Ke = obj.Ke * sparse(1:xNe, 1:xNe, obj.gea / 3);
                    % calculation of fi
                    obj.Fe = [zeros(3, xNe); repmat((obj.gea / 3), 3, 1)];
            end
            
            disp('************************************')
            disp('Calculation of Ke, Fe and C completed.')
            toc
        end
        
        function evalMe_TL3(obj)
            if obj.isMe_TL3_Evaluated, return; end
            tic, disp('-------------------------------------------------------');
            obj.mtcs.Me = repmat((1/6) * [1; 0.5; 1; 0.5; 0.5; 1], 1, obj.Ne) ...
                * sparse(1:obj.Ne, 1:obj.Ne, obj.gea);
            disp('Calculation of [Me] completed.');
            toc, disp('-------------------------------------------------------');
            % changing states
            obj.isMe_TL3_Evaluated = true;
        end
                
        %% Higher Order Elements
        function gd2elements(obj)

            % check the base element type
            if ~isequal(obj.etype, 'TL3')
                error('For generating higher order elements, the base element type must be TL3');
            end

            % check if already second order elements are generated
            if obj.isd2ElementsGenerated, return; end
            obj.ggmesh;
            
            % add new nodes
            Nn_old = obj.Nn;
            obj.nodes = [obj.nodes; (obj.nodes(obj.edges(:, 1), :) + obj.nodes(obj.edges(:, 2), :)) / 2];
            obj.cl = [obj.cl, abs(obj.elements(:, 1:3)) + Nn_old];
            obj.edges(:,end+1) = (1:size(obj.edges,1))' + Nn_old;

            % settig the new element type
            obj.etype = 'TL6';

            % change states
            obj.isd2ElementsGenerated = true;

        end
        
        function gd3elements(obj)

            % check the base element type
            if ~isequal(obj.etype, 'TL3')
                error('For generating higher order elements, the base element type must be TL3');
            end

            % check if already second order elements are generated
            if obj.isd3ElementsGenerated, return; end
            obj.ggmesh;
            
            % add new nodes
            Nn_old = obj.Nn;
            Nedges = size(obj.edges, 1);

            obj.nodes = [obj.nodes; ...
                obj.nodes(obj.edges(:,1),:) * 1/3 + obj.nodes(obj.edges(:,2),:) * 2/3; ...
                obj.nodes(obj.edges(:,1),:) * 2/3 + obj.nodes(obj.edges(:,2),:) * 1/3; ...
                (obj.nodes(obj.cl(:,1),:) + obj.nodes(obj.cl(:,2),:) + obj.nodes(obj.cl(:,3),:))/3];
            
            obj.cl = [obj.cl, abs(obj.elements(:, 1:3)) + Nn_old, abs(obj.elements(:, 1:3)) + Nn_old + Nedges, (1:obj.Ne)' + Nn_old + 2*Nedges];

            % settig element type
            obj.etype = 'TL10';
            % change states
            obj.isd3ElementsGenerated = true;
            
        end
        
        %% Topological Functions
        % setting needed data
        function setdata(obj)
            % first edge of each triangle
            e1 = obj.cl(:, [1, 2]);
            % second edge of each triangle
            e2 = obj.cl(:, [2, 3]);
            % third edge of each triangle
            e3 = obj.cl(:, [3, 1]);
            % sorting for lower index
            [e1, s1] = sort(e1, 2);
            [e2, s2] = sort(e2, 2);
            [e3, s3] = sort(e3, 2);
            % specefying chaNed edge index
            s1 = s1(:, 1) == 2;
            s2 = s2(:, 1) == 2;
            s3 = s3(:, 1) == 2;
            % unification of edges
            [obj.edges, ~, ic] = unique([e1; e2; e3], 'rows');
            % getting number of elements
            ne = obj.Ne;
            % getting index of edge corresponding to each elements
            e1 = ic(1:ne);
            e2 = ic(1 + ne:2 * ne);
            e3 = ic(1 + 2 * ne:3 * ne);
            % specefying boundary edges
            obj.bedges = sparse([e1, e2, e3], ones(3 * ne, 1), ones(3 * ne, 1));
            obj.bedges = full(obj.bedges == 1);
            % specefying trace direction
            e1(s1) = -e1(s1);
            e2(s2) = -e2(s2);
            e3(s3) = -e3(s3);
            % element matrix
            obj.elements(:, 1:3) = [e1, e2, e3];
            % edge element
            obj.edges = [obj.edges, zeros(size(obj.edges, 1), 6)];
            tmdbc_evalee(obj.edges, obj.elements);
            obj.edgeNamedSelections.('none') = find(obj.bedges);
        end
        
        function strefine(obj)
            mzNames = fieldnames(obj.mzs);
            
            for i = 1:numel(mzNames)
                obj.mzs.(mzNames{i}).strefine;
            end
            
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function adrefine(obj, maxLength, aspectRatio)
            obj.ggmesh;
            % triangles that nedded to be refined
            % refine index
            refineFlag = false(obj.Ne, 1);
            baseEdge = zeros(obj.Ne, 1);
            
            for i = 1:obj.Ne
                [refineFlag(i), baseEdge(i)] = needRefine(i);
            end
            
            % while loop until two mesh constraint be satisfy
            while any(refineFlag)
                ri = find(refineFlag);
                
                for i = 1:length(ri)
                    bisectTriangle(ri);
                end
                
            end
            
            function bisectTriangle(ti)
                eIndex = abs(obj.elements(ti, lIndex(ti)));
                eRow = obj.edges(eIndex, :);
                
                % if edge is a boundary edge
                if obj.bedges(eIndex)
                    % new node
                    ptmp = (obj.nodes(eRow(1), :) + obj.nodes(eRow(2), :)) / 2;
                    % modifying nodes
                    obj.nodes(end + 1, :) = ptmp;
                    nIndex = obj.Nn;
                    % modifying edges
                    obj.edges(end + 1, :) = obj.edges(eIndex, :);
                    obj.edges(eIndex, 2) = nIndex;
                    obj.edges(end, 1) = nIndex;
                    obj.edges(end + 1, :) = obj.edges(end);
                    % modifying cl
                    
                    tmp = obj.cl(ti, :);
                    
                    switch baseEdge(ti)
                        case 1
                            obj.cl(ti, :) = [tmp(1), nIndex, tmp(3)];
                            obj.elements(ti, 1:3) = 0;
                            obj.cl(end + 1, :) = [nIndex, tmp(2), tmp(3)];
                        case 2
                            obj.cl(ti, :) = [tmp(1), tmp(2), nIndex];
                            obj.cl(end + 1, :) = [tmp(1), nIndex, tmp(3)];
                        case 3
                            obj.cl(ti, :) = [tmp(1), tmp(2), nIndex];
                            obj.cl(end + 1, :) = [tmp(2), tmp(3), nIndex];
                    end
                    
                    e1 = eRow;
                    e1(2) = nIndex;
                    e2 = eRow;
                    e2(1) = nIndex;
                    e3 = eRow;
                    disp('salam');
                else
                    
                end
                
            end
            
            function [y1, y2] = needRefine(ti)
                % edge length
                el1 = norm(obj.nodes(obj.cl(ti, 1), :) - obj.nodes(obj.cl(ti, 2), :));
                el2 = norm(obj.nodes(obj.cl(ti, 2), :) - obj.nodes(obj.cl(ti, 3), :));
                el3 = norm(obj.nodes(obj.cl(ti, 3), :) - obj.nodes(obj.cl(ti, 1), :));
                [elMax, indexMax] = max([el1, el2, el3]);
                elMin = min([el1, el2, el3]);
                
                if (elMax / elMin > aspectRatio) || (elMax > maxLength)
                    y1 = true;
                else
                    y1 = false;
                end
                
                y2 = abs(obj.elements(ti, indexMax));
            end
            
        end
        
        %% Mesh Visiualization
        % show global mesh
        function varargout = showm(obj)

            obj.ggmesh;
            mzNames = fieldnames(obj.mzs);
            f = emdlab_r2d_mesh();
            ax = axes(f);
            f.Name = '[Global Mesh]';
            
            for i = 1:numel(mzNames)
                patch('Faces', obj.mzs.(mzNames{i}).cl, ...
                    'Vertices', obj.mzs.(mzNames{i}).nodes, 'FaceColor', ...
                    'c', 'EdgeColor', [0.2, 0.2, 0.2], ...
                    'FaceAlpha', 0.8, 'Parent', ax, 'ButtonDownFcn', @selectPatchCallbackGM);
            end

            index = obj.edges(:, 3) ~= obj.edges(:, 4);
            patch('Faces', obj.edges(index, [1, 2]), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2, 'parent', ax);
            
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
        
        % show free boundary
        function varargout = showfb(obj)

            obj.ggmesh;
            f = emdlab_r2d_mesh();
            ax = axes(f);
            f.Name = ['[Global Mesh Boundary Edges][', 'Nbe = ', num2str(sum(obj.bedges)), ']'];
            
            patch('Faces', obj.edges(obj.bedges, [1, 2]), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5, 'parent', ax);
            
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
        
        % show wire frame
        function varargout = showwf(obj)
            
            obj.ggmesh;
            f = emdlab_r2d_mesh();
            ax = axes(f);
            f.Name = '[Wire Frame of Mesh]';
            
            index = obj.edges(:, 3) ~= obj.edges(:, 4);
            patch('Faces', obj.edges(index, [1, 2]), 'Vertices', obj.nodes, ...
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
        
        % show mesh zones
        function varargout = showmzs(obj)
            
            mzNames = fieldnames(obj.mzs);
            f = emdlab_r2d_mesh();
            ax = axes(f);
            f.Name = ['[Number of Mesh Zones: ', num2str(numel(mzNames)), ']'];
            
            for i = 1:numel(mzNames)
                mzptr = obj.mzs.(mzNames{i});
                patch('Faces', mzptr.cl, ...
                    'Vertices', mzptr.nodes, 'FaceColor', ...
                    mzptr.color, 'linewidth', 0.05 ,'EdgeColor', [0, 0, 0], ...
                    'FaceAlpha', 1, 'Parent', ax);
            end
            
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
        
        % show geometry
        function varargout = showg(obj)

            obj.ggmesh;
            mzNames = fieldnames(obj.mzs);
            f = emdlab_r2d_mesh();
            ax = axes(f);
            f.Name = ['[Number of Zones: ', num2str(numel(mzNames)), ']'];
            
            for i = 1:numel(mzNames)
                mzptr = obj.mzs.(mzNames{i});
                p = patch('Faces', mzptr.cl, 'Vertices', mzptr.nodes, 'FaceColor', ...
                    mzptr.color, 'EdgeColor', 'none', ...
                    'FaceAlpha', 1, 'Parent', ax, 'ButtonDownFcn', @selectPatchCallback);
                p.UserData.color = mzptr.color;
            end
            
            index = obj.edges(:, 3) ~= obj.edges(:, 4);
            patch('Faces', obj.edges(index, [1, 2]), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 0.1, 'parent', ax, 'FaceAlpha', 0.95);

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
                       
        % show mesh degree
        function varargout = showmd(obj)

            obj.ggmesh;
            f = emdlab_r2d_mesh();
            ax = axes(f);
            f.Name = '[Global Mesh]';
            
            mzNames = fieldnames(obj.mzs);
            ecolor = zeros(obj.Ne, 3);
            
            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                ecolor(obj.ezi(:, i), 1) = mzptr.color(1);
                ecolor(obj.ezi(:, i), 2) = mzptr.color(2);
                ecolor(obj.ezi(:, i), 3) = mzptr.color(3);
            end
            
            % point index
            if isequal(obj.etype, 'TL3')
                pIndex = 1:3;
            elseif isequal(obj.etype, 'TL6')
                pIndex = 1:6;
            elseif isequal(obj.etype, 'TL10')
                pIndex = 1:10;
            end
            
            patch('Faces', obj.cl(:, [1,2,3]), 'Vertices', obj.nodes, ...
                'FaceColor', 'flat', 'FaceVertexCData', ecolor, ...
                'FaceAlpha', 0.5, 'EdgeColor', 'k', 'parent', ax);

            patch('Faces', obj.cl(:, pIndex), 'Vertices', obj.nodes, ...
                'FaceColor', 'none', 'EdgeColor', 'none', 'parent', ax, 'Marker', 'o', 'MarkerFaceColor', 'k');
            
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
        
        % show nodes on global mesh
        function varargout = showNodes(obj, varargin)

            [f,ax] = obj.showm;
            for i = 1:numel(varargin)
                patch('Faces', varargin{i}, 'Vertices', obj.nodes, ...
                    'EdgeColor', 'none', 'parent', ax, 'Marker', 'o', 'MarkerFaceColor', 'r');
            end

            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end
        
        function varargout = shownes(obj)
            f = GraphicWindow;
            f.Name = '[Named edges of mesh]';
            h = guihandles(f);
            nes = fieldnames(obj.edgeNamedSelections);
            
            for i = 1:numel(nes)
                tmp = rand(1, 3);
                patch('Faces', obj.edges(obj.edgeNamedSelections.(nes{i}), 1:2), 'Vertices', ...
                    obj.nodes, 'FaceColor', ...
                    tmp, 'EdgeColor', tmp, 'linewidth', 1.5, ...
                    'FaceAlpha', 1, 'parent', h.va);
            end
            
            legend(h.va, nes);
            set(f, 'HandleVisibility', 'off', 'Visible', 'on');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end
        
        function showne(obj, neName)
            neName = obj.checkEdgeNamedSelectionExistence(neName);
            ah = setFigure('TMDBC: Named Edge');
            patch('Faces', obj.edges(obj.edgeNamedSelections.(neName), 1:2), 'Vertices', ...
                obj.nodes, 'FaceColor', ...
                'b', 'EdgeColor', 'b', 'linewidth', 1.5, ...
                'FaceAlpha', 1, 'parent', ah);
            set(gcf, 'HandleVisibility', 'off', 'Visible', 'on');
        end
        
        function varargout = showvf(obj, Fx, Fy)
            f = GraphicWindow;
            f.Name = 'Field plot on center of elements';
            h = guihandles(f);
            h.va.NextPlot = 'add';
            c = obj.getCenterOfElements;
            color = zeros(obj.Ne, 3);
            mzNames = fieldnames(obj.mzs);
            
            for i = 1:obj.Nmzs
                mzptr = obj.mzs.(mzNames{i});
                color(obj.ezi(:, mzptr.zi), :) = repmat(mzptr.color, mzptr.Ne, 1);
            end
            
            patch(h.va, 'faces', obj.cl(:, 1:3), 'vertices', obj.nodes, ...
                'facecolor', 'flat', 'FaceVertexCData', color, 'EdgeColor', 'w', 'facealpha', 0.5);
            quiver(h.va, c(:, 1), c(:, 2), Fx, Fy, 'color', 'k');
            axis(h.va, 'off', 'equal');
            set(f, 'Visible', 'on');
            
            if nargout == 1
                varargout{1} = f;
            elseif nargout > 1
                error('Too many output argument.');
            end
            
        end
        
        %% Tools: copy and transform
        % copy and transform
        function copyMirrorMeshZone(obj, nmzName, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getMirror(varargin{:});
        end
        
        function cmmz(varargin)
            copyMirrorMeshZone(varargin{:});
        end
        
        function copyRotateMeshZone(obj, nmzName, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getRotate(varargin{:});
        end
        
        function crmz(varargin)
            copyRotateMeshZone(varargin{:});
        end
        
        function copyShiftMeshZone(obj, nmzName, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            obj.mzs.(nmzName) = obj.mzs.(mzName).getShift(varargin{:});
        end
        
        function cshmz(varargin)
            copyShiftMeshZone(varargin{:});
        end
        
        % only transform
        function mirrorMeshZone(obj, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            obj.mzs.(mzName).mirror(varargin{:});
        end
        
        function mmz(varargin)
            mirrorMeshZone(varargin{:});
        end
        
        function rotateMeshZone(obj, mzName, varargin)
            mzName = obj.checkMeshZoneExistence(mzName);
            obj.mzs.(mzName).rotate(varargin{:});
        end
        
        function rmz(varargin)
            rotateMeshZone(varargin{:});
        end
        
        function shiftMeshZone(obj, mzName, varargin)
            mzName = rmspaces(mzName);
            obj.checkMeshZoneExistence(mzName)
            obj.mzs.(mzName).shift(varargin{:});
        end
        
        function shmz(varargin)
            shiftMeshZone(varargin{:});
        end
        
        %% Tools: some operations on mesh zones
        function joinMeshZones(obj, nmzName, varargin)
            xNmzs = numel(varargin);
            
            if xNmzs < 2
                error('Minimum number mzs must be 2.');
            end
            
            for i = 1:xNmzs
                varargin{i} = obj.checkMeshZoneExistence(varargin{i});
            end
            
            nmzName = obj.checkMeshZoneNonExistence(nmzName);
            Nn_tmp = zeros(1, xNmzs);
            Ne_tmp = zeros(1, xNmzs);
            
            for i = 1:numel(varargin)
                Nn_tmp(i) = obj.mzs.(varargin{i}).Nn;
                Ne_tmp(i) = obj.mzs.(varargin{i}).Ne;
            end
            
            n_nmz = zeros(sum(Nn_tmp), 2);
            e_nmz = zeros(sum(Ne_tmp), 3);
            n_tmp = 0;
            e_tmp = 0;
            
            for i = 1:numel(varargin)
                n_nmz(1 + n_tmp:n_tmp + Nn_tmp(i), :) = obj.mzs.(varargin{i}).nodes;
                e_nmz(1 + e_tmp:e_tmp + Ne_tmp(i), :) = obj.mzs.(varargin{i}).cl + n_tmp;
                n_tmp = n_tmp + Nn_tmp(i);
                e_tmp = e_tmp + Ne_tmp(i);
            end
            
            % jointing mzs
            [n_nmz, ~, ic] = uniquetol(n_nmz, obj.gleps, 'ByRows', true);
            e_nmz = ic(e_nmz);
            % adding new mz
            obj.mzs.(nmzName) = emdlab_m2d_tmz(e_nmz, n_nmz);
            obj.mzs.(nmzName).material = obj.mzs.(varargin{1}).material;
            obj.mzs.(nmzName).color = obj.mzs.(varargin{1}).color;
            % removing old mzs
            for i = 1:numel(varargin)
                obj.mzs = rmfield(obj.mzs, varargin{i});
            end
            
        end
        
        function jmzs(varargin)
            joinMeshZones(varargin{:});
        end
        
        function getQuality(obj)
            obj.ggmesh;
            % edges length
            el = sqrt(sum((obj.nodes(obj.edges(:, 1), :) - ...
                obj.nodes(obj.edges(:, 2), :)).^2, 2));
            b1 = el(abs(obj.elements(:, 1)));
            b2 = el(abs(obj.elements(:, 2)));
            b3 = el(abs(obj.elements(:, 3)));
            % mesh quality
            y = ((b1 + b2 - b3) .* (b1 - b2 + b3) .* (-b1 + b2 + b3)) ./ (b1 .* b2 .* b3);
            fprintf('Average Quality = %f\n', mean(y));
            fprintf('Minimum Quality = %f\n', min(y));
        end
        
        function gq(varargin)
            getQuality(varargin{:});
        end
        
        function ttmdbc = getExtrude(obj, z, skewAngle)
            
            ttmdbc = emdlab_m3d_ttmdb;
            mzNames = fieldnames(obj.mzs);
            
            if iscolumn(z)
                z = z';
            end
            
            Nz = length(z);
            
            if nargin < 3
                
                for i = 1:numel(mzNames)
                    mzptr = obj.mzs.(mzNames{i});
                    ztmp = repmat(z, mzptr.Nn, 1);
                    ttmdbc.addmz(mzNames{i}, emdlab_m3d_ttmz(tmzpc_getExtrude(...
                        mzptr.cl, obj.elements(obj.ezi(:, mzptr.zi), 1:3), ...
                        mzptr.Nn, Nz - 1), [repmat(mzptr.nodes, Nz, 1), ztmp(:)]));
                    ttmdbc.mzs.(mzNames{i}).material = mzptr.material;
                end
                
            else
                stepAngle = skewAngle / (Nz - 1);
                p = zeros(obj.Nn * Nz, 3);
                
                for i = 1:Nz
                    p((i - 1) * obj.Nn + 1:i * obj.Nn, 1:3) = ...
                        ext_protate3z([obj.nodes, repmat(z(i), obj.Nn, 1)], (i - 1) * stepAngle);
                end
                
                ttmz = emdlab_m3d_ttmz(tmzpc_getExtrude(obj.cl, obj.elements, ...
                    obj.Nn, Nz - 1), p);
            end
            
        end
        
        function getMeshZoneExtrude(obj, ttmptr, mzName, z, skewAngle)
            mzName = obj.checkMeshZoneExistence(mzName);
            
            if iscolumn(z)
                z = z';
            end
            
            Nz = length(z);
            mzptr = obj.mzs.(mzName);
            ztmp = repmat(z, mzptr.Nn, 1);
            
            if nargin < 5
                ttmptr.addmz(mzName, emdlab_m3d_ttmz(tmzpc_getExtrude(...
                    mzptr.cl, obj.elements(obj.ezi(:, mzptr.zi), 1:3), ...
                    mzptr.Nn, Nz - 1), [repmat(mzptr.nodes, Nz, 1), ztmp(:)]));
            else
                stepAngle = skewAngle * (pi / 180) / (Nz - 1);
                p = zeros(mzptr.Nn * Nz, 3);
                
                for i = 1:Nz
                    p((i - 1) * mzptr.Nn + 1:i * mzptr.Nn, 1:3) = ...
                        ext_protate3z([mzptr.nodes, repmat(z(i), mzptr.Nn, 1)], (i - 1) * stepAngle);
                end
                
                ttmptr.addmz(mzName, TTMZPC(tmzpc_getExtrude(...
                    mzptr.cl, obj.elements(obj.ezi(:, mzptr.zi), 1:3), ...
                    mzptr.Nn, Nz - 1), p));
            end
            
            ttmptr.mzs.(mzName).material = mzptr.material;
        end
        
        %% Index Operations
        % free boundary edge and node index
        function y = getfbe(obj)
            obj.ggmesh;
            y = find(obj.bedges);
        end
        
        function y = getfbn(obj)
            y = obj.getfbe;
            if obj.etype == 'TL3'
            y = obj.edges(y, [1,2]);
            elseif obj.etype == 'TL6'
            y = obj.edges(y, [1,2,end]);
            end
            y = unique(y(:));
        end
        
        % free boundary edge index on line
        function y = getfbeiol_p0u(obj, varargin)
            obj.ggmesh;
            y = tmdbc_getfbeiol(obj.edges, obj.nodes, obj.edgeNamedSelections.none, varargin{:});
            y = y(y ~= 0);
        end
        
        function y = getfbeiol_p0p1(obj, p0, p1)
            y = obj.getfbeiol_p0u(p0, p1 - p0);
        end
        
        function y = getfbeiol_p0x(obj, p0)
            y = obj.getfbeiol_p0u(p0, [1, 0]);
        end
        
        function y = getfbeiol_p0y(obj, p0)
            y = obj.getfbeiol_p0u(p0, [0, 1]);
        end
        
        function y = getfbniol_p0u(obj, varargin)
            % get|freebounday|node|index|on|line
            y = getfbeiol_p0u(obj, varargin{:});
            if obj.etype == 'TL3'
            y = obj.edges(y, [1,2]);
            elseif obj.etype == 'TL6'
            y = obj.edges(y, [1,2,end]);
            end
            y = unique(y(:));
        end
        
        function y = getfbniol_p0p1(obj, p0, p1)
            y = obj.getfbniol_p0u(p0, p1 - p0);
        end
        
        function y = getfbniol_p0x(obj, p0)
            y = obj.getfbniol_p0u(p0, [1, 0]);
        end
        
        function y = getfbniol_p0y(obj, p0)
            y = obj.getfbniol_p0u(p0, [0, 1]);
        end
        
        % respect to circle
        function y = getfbeioc(obj, varargin)
            obj.ggmesh;
            y = tmdbc_getfbeioc(obj.edges, obj.nodes, obj.edgeNamedSelections.none, varargin{:});
            y = y(y ~= 0);
        end
        
        function y = getfbeiic(obj, varargin)
            obj.ggmesh;
            y = tmdbc_getfbeiic(obj.edges, obj.nodes, obj.edgeNamedSelections.none, varargin{:});
            y = y(y ~= 0);
        end
        
        function y = getfbeioutc(obj, varargin)
            obj.ggmesh;
            y = tmdbc_getfbeioutc(obj.edges, obj.nodes, obj.edgeNamedSelections.none, varargin{:});
            y = y(y ~= 0);
        end
        
        function y = getfbnioc(obj, varargin)
            y = getfbeioc(obj, varargin{:});
            if obj.etype == 'TL3'
            y = obj.edges(y, [1,2]);
            elseif obj.etype == 'TL6'
            y = obj.edges(y, [1,2,end]);
            end
            y = unique(y(:));
        end
        
        function y = getfbniic(obj, varargin)
            y = getfbeiic(obj, varargin{:});
            y = obj.edges(y, 1:2);
            y = unique(y(:));
        end
        
        function y = getfbnioutc(obj, varargin)
            y = getfbeioutc(obj, varargin{:});
            y = obj.edges(y, 1:2);
            y = unique(y(:));
        end
        
        % node index
        function y = getnIndexOnPoint(obj, x, y, tol)
            
            if nargin < 4
                tol = obj.gleps;
            end
            
            y = [obj.nodes(:, 1) - x, obj.nodes(:, 2) - y];
            y = sqrt(sum(y.^2, 2));
            y = find(y < tol);
        end

        function y = getnIndexOnCircle(obj, c, r, tol)
            
            if nargin < 4
                tol = obj.gleps;
            end
            
            y = [obj.nodes(:, 1) - c(1), obj.nodes(:, 2) - c(2)];
            y = sqrt(sum(y.^2, 2));
            y = find(abs(y - r) < tol);
        end
        
        function y = getnIndexOnLine(obj, p1, p2, tol)
            
            if nargin < 4
                tol = obj.gleps;
            end
            
            u = (p2 - p1) / norm(p2 - p1);
            pp1 = [obj.nodes(:, 1) - p1(1), obj.nodes(:, 2) - p1(2)];
            alpha = (pp1 * u');
            y = sqrt(sum((pp1 - alpha * u).^2, 2));
            y = find(y < tol);
        end
        
        function y = getnIndexOnRay(obj, p1, p2)
            u = (p2 - p1) / norm(p2 - p1);
            pp1 = [obj.nodes(:, 1) - p1(1), obj.nodes(:, 2) - p1(2)];
            alpha = (pp1 * u');
            y = sqrt(sum((pp1 - alpha * u).^2, 2));
            y = find(y < obj.gleps);
            index = alpha(y) <- obj.gleps;
            y = y(~index);
        end
        
        function y = getnIndexOnEdges(obj, eList)
            
            switch obj.etype
                case 'TL3'
                    y = obj.edges(eList, 1:2);
                    y = unique(y(:));
                case 'TL6'
                    y = obj.edges(eList, 1:2);
                    y = unique(y(:));
                    y = [y; eList + obj.Nn - size(obj.edges, 1)];
            end
            
        end
        
        % edge index
        function y = geteIndexOnCircle(obj, varargin)
            y = obj.getnIndexOnCircle(varargin{:});
            y = bitand(ismember(obj.edges(:, 1), y), ismember(obj.edges(:, 2), y));
            y = find(y);
        end
        
        function y = geteIndexOnLine(obj, varargin)
            y = obj.getnIndexOnLine(varargin{:});
            y = bitand(ismember(obj.edges(:, 1), y), ismember(obj.edges(:, 2), y));
            y = find(y);
        end
        
        function y = getCenterOfElements(obj)
            obj.ggmesh;
            y = (obj.nodes(obj.cl(:, 1), :) + obj.nodes(obj.cl(:, 2), :) + obj.nodes(obj.cl(:, 3), :)) / 3;
        end
        
        % periodic nodes
        function [km, ks] = splitRotate(obj, k, varargin)
            Nk = length(k);
            
            if rem(Nk, 2) ~= 0
                error('number of input indices must be even.');
            end
            
            polarAngle = atan_02pi(obj.nodes(k, :));
            [~, index] = sort(polarAngle);
            k = k(index);
            pm = obj.nodes(k(1:Nk / 2), :);
            ps = obj.nodes(k(Nk / 2 + 1:Nk), :);
            km = k(1:Nk / 2);
            ks = k(Nk / 2 + 1:Nk);
            [Flag, index] = ismembertol(ps, pm, obj.gleps, 'ByRow', true);
            
            if any(~Flag)
                error('These set of points does not form a set of peridic points.');
            end
            
            ks = ks(index);
        end
        
        function [km, ks] = splitPeriodic(obj, k, varargin)
            Nk = length(k);
            
            if rem(Nk, 2) ~= 0
                error('number of input indices must be even.');
            end
            
            tp = obj.nodes(k, :);
            
            temp = 1;
            
            while ~isempty(tp)
                sp = tp(1, :);
                tp = tp(2:end, :);
                km(temp) = k(1);
                k = k(2:end);
                sp = ext_protate2(sp, varargin{1});
                index = find(sqrt(sum([tp(:, 1) - sp(1), tp(:, 2) - sp(2)].^2, 2)) < obj.gleps);
                
                if index
                    ks(temp) = k(index);
                    tp = tp([1:index - 1, index + 1:end], :);
                    k = k([1:index - 1, index + 1:end]);
                    temp = temp + 1;
                    continue
                end
                
                sp = ext_protate2(sp, -2 * varargin{1});
                index = find(sqrt(sum([tp(:, 1) - sp(1), tp(:, 2) - sp(2)].^2, 2)) < obj.gleps);
                
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
        
        function [km, ks] = splitShift(obj, k, varargin)
            % Number of points
            Nk = length(k);
            
            if rem(Nk, 2) ~= 0
                error('Number of points indices must be an even number.');
            end
            
            % total points
            tp = obj.nodes(k, :);
            
            temp = 1;
            
            while ~isempty(tp)
                sp = tp(1, :);
                tp = tp(2:end, :);
                km(temp) = k(1);
                k = k(2:end);
                sp = ext_pshift2(sp, varargin{1});
                index = find(sqrt(sum([tp(:, 1) - sp(1), tp(:, 2) - sp(2)].^2, 2)) < obj.gleps);
                
                if index
                    ks(temp) = k(index);
                    tp = tp([1:index - 1, index + 1:end], :);
                    k = k([1:index - 1, index + 1:end]);
                    temp = temp + 1;
                    continue
                end
                
                sp = ext_pshift2(sp, -2 * varargin{1});
                index = find(sqrt(sum([tp(:, 1) - sp(1), tp(:, 2) - sp(2)].^2, 2)) < obj.gleps);
                
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
        
        %% Named Selections
        % edge
        function name = checkEdgeNamedSelectionExistence(obj, name)
            name = rmspaces(name);
            
            if ~isfield(obj.edgeNamedSelections, name)
                error('Specified edge named selection does not exist.');
            end
            
        end
        
        function name = checkEdgeNamedSelectionNonExistence(obj, name)
            name = rmspaces(name);
            
            if isfield(obj.edgeNamedSelections, name)
                error('Specified edge named selection already exist.');
            end
            
        end
        
        function addEdgeNamedSelection(obj, name, indices)
            name = obj.checkEdgeNamedSelectionNonExistence(name);
            
            if isempty(indices)
                error('indices matreix is empty.');
            end
            
            obj.edgeNamedSelections.(name) = indices;
            obj.edgeNamedSelections.('none') = setdiff(...
                obj.edgeNamedSelections.('none'), ...
                obj.edgeNamedSelections.(name));
        end
        
        %% Import and Export
        function read_g2d_txt(obj, Dir, MG)
            
            if nargin < 3
                MG = 'MM';
            else
                MG = upper(strrep(MG, ' ', ''));
            end
            
            f = fopen(Dir, 'r');
            clc
            tic
            
            while true
                str = fgetl(f);
                % check for end of file
                if isnumeric(str), break; end
                str = strsplit(str, ' ');
                
                if strcmpi(str{1}, 'BeginFace:')
                    mzName = str{2};
                    str = strsplit(fgetl(f), ' ');
                    Nloop = str2double(str{1});
                    f1 = cell(1, Nloop);
                    f2 = cell(1, Nloop);
                    v1 = cell(1, Nloop);
                    v2 = cell(1, Nloop);
                    NoldNodes = 0;
                    
                    for i = 1:Nloop
                        str = strsplit(fgetl(f), ' ');
                        
                        if strcmpi(str{1}, 'BeginLoop:')
                            % reading connectivity
                            str = fgetl(f);
                            Nf = str2double(str);
                            f1{i} = zeros(1, Nf);
                            f2{i} = zeros(1, Nf);
                            
                            for j = 1:Nf
                                str = fgetl(f);
                                str = strsplit(str);
                                f1{i}(j) = str2double(str{1});
                                f2{i}(j) = str2double(str{2});
                            end
                            
                            if i > 1
                                NoldNodes = NoldNodes + length(v1{i - 1});
                                f1{i} = f1{i} + NoldNodes;
                                f2{i} = f2{i} + NoldNodes;
                            end
                            
                            % reading vertices
                            str = fgetl(f);
                            Nv = str2double(str);
                            v1{i} = zeros(1, Nv);
                            v2{i} = zeros(1, Nv);
                            
                            for j = 1:Nv
                                str = fgetl(f);
                                str = strsplit(str);
                                v1{i}(j) = str2double(str{1});
                                v2{i}(j) = str2double(str{2});
                            end
                            
                            str = fgetl(f);
                            str = strsplit(str);
                            
                            if ~strcmpi(str{1}, 'EndLoop.')
                                error('Wrong File Type.');
                            end
                            
                        end
                        
                    end
                    
                    str = fgetl(f);
                    str = strsplit(str);
                    
                    if ~strcmpi(str{1}, 'EndFace.')
                        error('Wrong File Type.');
                    end
                    
                end
                
                % add to mesh database
                con = [cell2mat(f1); cell2mat(f2)]';
                ver = [cell2mat(v1); cell2mat(v2)]';
                [ver, ~, ic] = uniquetol(ver, 1e-6, 'ByRows', true);
                con = ic(con);
                
                switch upper(MG)
                    case 'MM'
                        obj.addmz(mzName, MinimalMesh(con, ver));
                    case 'MG0'
                        obj.addmz(mzName, MeshGenerator0(con, ver));
                    case 'MG1'
                        obj.addmz(mzName, MeshGenerator1(con, ver));
                    case 'MG3'
                        obj.addmz(mzName, MeshGenerator3(con, ver));
                    otherwise
                        error('Specified mesh generator does not exist.');
                end
                
            end
            
            fclose(f);
            toc
        end
        
        function read_g2d_bin(obj, Dir, MG, varargin)
            
            if nargin < 3
                MG = 'MM';
            end
            
            f = fopen(Dir, 'r');
            clc
            tic
            
            while true
                str = fread(f, 80, '*char');
                
                if isempty(str)
                    break;
                end
                
                str = strsplit(strtrim(str'));
                
                if strcmp(str{1}, 'BeginFace:')
                    mzname = str{2};
                    Nloop = fread(f, 1, 'uint32');
                    f1 = cell(1, Nloop);
                    f2 = cell(1, Nloop);
                    v1 = cell(1, Nloop);
                    v2 = cell(1, Nloop);
                    NoldNodes = 0;
                    
                    for i = 1:Nloop
                        str = fread(f, 80, '*char');
                        str = strsplit(strtrim(str'));
                        
                        if strcmp(str{1}, 'BeginLoop:')
                            % reading connectivity
                            Nf = fread(f, 1, 'uint32');
                            f1{i} = fread(f, Nf, 'uint32')';
                            f2{i} = fread(f, Nf, 'uint32')';
                            
                            if i > 1
                                NoldNodes = NoldNodes + length(v1{i - 1});
                                f1{i} = f1{i} + NoldNodes;
                                f2{i} = f2{i} + NoldNodes;
                            end
                            
                            % reading vertices
                            Nv = fread(f, 1, 'uint32');
                            v1{i} = fread(f, Nv, 'double')';
                            v2{i} = fread(f, Nv, 'double')';
                            str = fread(f, 80, '*char');
                            str = strsplit(strtrim(str'));
                            
                            if ~strcmp(str{1}, 'EndLoop')
                                error('Wrong File Type.');
                            end
                            
                        end
                        
                    end
                    
                    str = fread(f, 80, '*char');
                    str = strsplit(strtrim(str'));
                    
                    if ~strcmp(str{1}, 'EndFace')
                        error('Wrong File Type.');
                    end
                    
                end
                
                % add to mesh database
                con = [cell2mat(f1); cell2mat(f2)]';
                ver = [cell2mat(v1); cell2mat(v2)]';
                [ver, ~, ic] = uniquetol(ver, 1e-6, 'ByRows', true);
                con = ic(con);
                
                switch upper(MG)
                    case 'MM'
                        obj.addmz(mzname, emdlab_m2d_mm(con, ver, varargin{:}));
                    case 'MM1'
                        obj.addmz(mzname, emdlab_m2d_mm(con, ver, varargin{:}));
                    case 'MG0'
                        obj.addmz(mzname, emdlab_m2d_mg0(con, ver, varargin{:}));
                    case 'MG1'
                        obj.addmz(mzname, emdlab_m2d_mg1(con, ver, varargin{:}));
                    case 'MG3'
                        obj.addmz(mzname, emdlab_m2d_mg3(con, ver, varargin{:}));
                    otherwise
                        error('Invalid mesh generator name.');
                end
                
            end
            
            fclose(f);
            toc
        end
        
        function read_stlf_bin(obj, Dir)
            
            if ~isdir(Dir)
                error('Directory was not found.');
            end
            
            stlFile = ls(Dir);
            stlFile = stlFile(3:end, :);
            
            for i = 1:size(stlFile, 1)
                obj.read_stl_bin([Dir, '\', stlFile(i, :)]);
            end
            
        end
        
        function read_stl_bin(obj, Dir)
            % openning STL file
            f = fopen(Dir, 'r');
            
            if f < 0
                error('Can not open file.');
            end
            
            % zone name
            mzname = fread(f, 80, '*char');
            mzname = strrep(mzname', ' ', '');
            % number of elements
            Nt = fread(f, 1, 'uint32');
            % trianle points
            p = zeros(9, Nt);
            % loop over elements
            for i = 1:Nt
                % reading normal vector
                fread(f, 3, 'single');
                % reading vertices
                p(:, i) = fread(f, 9, 'single');
                % reading color
                fread(f, 1, 'uint16');
            end
            
            % closing STL file
            fclose(f);
            % generation of new mesh zone
            % nodes
            p = reshape(p, 3, []);
            p = p';
            % elements
            e = 1:3 * Nt;
            e = reshape(e, 3, []);
            e = e';
            % unification
            [p, ~, index] = uniquetol(p, 1e-6, 'ByRows', true);
            % reindex
            e = index(e);
            % adding mesh zone
            obj.addmz(mzname, emdlab_m2d_tmz(e, p(:, 1:2)));
        end
        
        function write_stl_bin(obj)
            % openning STL file
            if exist('STL_Files', 'file')
                rmdir('STL_Files', 's')
            end
            
            mkdir('STL_Files')
            mznames = fieldnames(obj.mzs);
            
            for i = 1:numel(mznames)
                % zone name
                mzname = mznames{i};
                f = fopen(['STL_Files\', mzname, '.STL'], 'w');
                % writing zone name
                fwrite(f, [mzname, repmat(' ', 1, 80 - length(mzname))], '*char');
                % number of elements
                fwrite(f, obj.mzs.(mzname).Ne, 'uint32');
                % loop over elements
                for j = 1:obj.mzs.(mzname).Ne
                    % writing normal vector
                    fwrite(f, [0, 0, 1], 'single');
                    % writing vertices
                    ti = obj.mzs.(mzname).cl(j, :);
                    fwrite(f, [obj.mzs.(mzname).nodes(ti(1), :), 0]', 'single');
                    fwrite(f, [obj.mzs.(mzname).nodes(ti(2), :), 0]', 'single');
                    fwrite(f, [obj.mzs.(mzname).nodes(ti(3), :), 0]', 'single');
                    % writing color
                    fwrite(f, 1, 'uint16');
                end
                
                % closing STL file
                fclose(f);
            end
            
        end
        
        function write_m2d_bin(obj, Dir)
            % openning STL foler
            if nargin < 2
                
                if exist('m2d_Files', 'file')
                    rmdir('m2d_Files', 's');
                end
                
                mkdir('m2d_Files');
                Dir = 'm2d_Files';
            end
            
            % loop over mesh zones
            mznames = fieldnames(obj.mzs);
            
            for i = 1:numel(mznames)
                % zone name
                mzname = mznames{i};
                f = fopen([Dir, '\', mzname, '.m2d'], 'w');
                % writing number of nodes and elements
                fwrite(f, uint32(obj.mzs.(mzname).Nn), 'uint32');
                fwrite(f, uint32(obj.mzs.(mzname).Ne), 'uint32');
                % writing nodes
                fwrite(f, obj.mzs.(mzname).nodes(:), 'double');
                % writing elements
                fwrite(f, obj.mzs.(mzname).cl(:), 'double');
                % closing m2d file
                fclose(f);
            end
            
        end
        
        function read_m2d_bin(obj, Dir)
            % openning m2d file
            f = fopen(Dir, 'r');
            Nnodes = fread(f, 1, 'uint32');
            Nelements = fread(f, 1, 'uint32');
            p = fread(f, 2 * Nnodes, 'double');
            e = fread(f, 3 * Nelements, 'double');
            fclose(f);
            p = reshape(p, [], 2);
            e = reshape(e, [], 3);
            % adding mesh zone
            mzname = strsplit(Dir, '\');
            mzname = strsplit(mzname{end}, '.');
            obj.addmz(mzname{1}, emdlab_m2d_tmz(e, p));
        end
        
        function read_m2df_bin(obj, Dir)
            
            if ~isdir(Dir)
                error('Directory was not found.');
            end
            
            stlFile = ls(Dir);
            stlFile = stlFile(3:end, :);
            
            for i = 1:size(stlFile, 1)
                obj.read_m2d_bin([Dir, '\', stlFile(i, :)]);
            end
            
        end
        
        function read_msh_bin(obj, Dir)
            % openning msh file
            f = fopen(Dir, 'r');
            
            if f < 0
                error('Can not open file.');
            end
            
            % allocating space for mesh zone elements
            e = cell(1, 20);
            
            for i = 1:numel(e)
                e{i} = zeros(3, []);
            end
            
            % loop for read
            while true
                str = fgetl(f);
                
                if strcmp(str, '$Nodes')
                    Nnodes = fgetl(f);
                    Nnodes = str2double(Nnodes);
                    p = zeros(3, Nnodes);
                    
                    for i = 1:Nnodes
                        fread(f, 1, 'int32');
                        p(:, i) = fread(f, 3, 'double');
                    end
                    
                end
                
                if strcmp(str, '$Elements')
                    Nnodes = fgetl(f);
                    Nnodes = str2double(Nnodes);
                    
                    for i = 1:Nnodes
                        eIndex = fread(f, 1, 'int32');
                        
                        if eIndex == 15
                            Nelements = fread(f, 1, 'int32');
                            Nt = fread(f, 1, 'int32');
                            
                            for j = 1:Nelements
                                fread(f, 1, 'int32');
                                fread(f, Nt, 'int32');
                                fread(f, 1, 'int32');
                            end
                            
                        elseif eIndex == 1
                            Nelements = fread(f, 1, 'int32');
                            Nt = fread(f, 1, 'int32');
                            
                            for j = 1:Nelements
                                fread(f, 1, 'int32');
                                fread(f, Nt, 'int32');
                                fread(f, 2, 'int32');
                            end
                            
                        elseif eIndex == 2
                            Nelements = fread(f, 1, 'int32');
                            Nt = fread(f, 1, 'int32');
                            
                            for j = 1:Nelements
                                fread(f, 1, 'int32');
                                fread(f, Nt - 1, 'int32');
                                zIndex = fread(f, 1, 'int32');
                                e{zIndex}(:, end + 1) = fread(f, 3, 'int32');
                            end
                            
                        end
                        
                    end
                    
                end
                
                if strcmp(str, '$EndElements'), break; end
            end
            
            fclose(f);
            % adding mesh zone
            e = e(1:zIndex);
            p = p(1:2, :)';
            
            for i = 1:zIndex
                xtmp = p(:, 1);
                ytmp = p(:, 2);
                xtmp = xtmp(e{i});
                ytmp = ytmp(e{i});
                [ptmp, ~, index] = uniquetol([xtmp(:), ytmp(:)], 1e-6, 'ByRows', true);
                etmp = 1:3 * size(e{i}, 2);
                etmp = etmp(index);
                etmp = reshape(etmp, 3, [])';
                obj.addmz(emdlab_m2d_tmz(etmp, ptmp));
            end
            
        end
        
        %% Auxiliary tools
        function makeFalse_isGlobalMeshGenerated(obj)

            obj.isGlobalMeshGenerated = false;
            obj.isd2ElementsGenerated = false;
            obj.isd3ElementsGenerated = false;
            obj.isJITEvaluated = false;
            obj.isKeMeFe_TL3_Evaluated = false;
            obj.isKeMeFe_TL6_Evaluated = false;
            obj.isKeFe_TL3_Evaluated = false;
            obj.isKeFe_TL6_Evaluated = false;
            obj.isKexy4Fe_TL3_Evaluated = false;
            obj.isKeFe_TL6_cte_Evaluated = false;
            obj.isKe_TL3_Evaluated = false;
            obj.isFe_TL3_Evaluated = false;
            obj.isMe_TL3_Evaluated = false;
            obj.isKexy1_TL3_Evaluated = false;
            obj.isKexy4_TL3_Evaluated = false;
            obj.isKexy1_TL6_Evaluated = false;
            obj.isKerz1_TL3_Evaluated = false;

        end
        
        function aux_crjmzs(obj, mzname, Nt, Ns)
            % copy|rotate|join|mesh|zones
            if nargin < 4, Ns = Nt; end
            
            tau_t = 2 * pi / Nt;
            tmp = obj.getDefaultMeshZoneName;
            
            obj.changeMeshZoneName(mzname, [tmp, '1']);
            
            for i = 1:Ns-1
                obj.crmz([tmp, num2str(i + 1)], [tmp, num2str(i)], tau_t);
            end
            
            tmp = getlist(tmp, 1:Ns);
            obj.jmzs(mzname, tmp{:});
            % change states
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function aux_cmrjmzs(obj, mzname, Nt, Ns)
            % copy|mirror|rotate|join
            if nargin < 4, Ns = Nt; end
            tau_t = 2 * pi / Nt;
            % mirror: slot axis
            mirrorAxis = [cos(tau_t / 2), sin(tau_t / 2)];
            tmp = obj.getDefaultMeshZoneName;
            obj.changeMeshZoneName(mzname, [tmp, '1']);
            obj.cmmz([tmp, '2'], [tmp, '1'], mirrorAxis);
            
            for i = 1:2:2 * (Ns - 1)
                obj.crmz([tmp, num2str(i + 2)], [tmp, num2str(i)], tau_t);
                obj.crmz([tmp, num2str(i + 3)], [tmp, num2str(i + 1)], tau_t);
            end
            
            tmp = getlist(tmp, 1:2 * Ns);
            obj.jmzs(mzname, tmp{:});
            % change states
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function aux_cshjmzs(obj, mzname, Nc, vec)
            % copy|shift|join|mesh|zones
            if Nc < 2, error('Number of copies must be grater than 2.'); end
            tmp = obj.getDefaultMeshZoneName;
            obj.changeMeshZoneName(mzname, [tmp, '1']);
            
            for i = 2:Nc
                obj.cshmz([tmp, num2str(i)], [tmp, num2str(i - 1)], vec);
            end
            
            tmp = getlist(tmp, 1:Nc);
            obj.jmzs(mzname, tmp{:});
            % change states
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function aux_cmrjmzx(obj, mzname, Nt, Ns)
            % copy|mirror|rotate|join
            if nargin < 4, Ns = Nt; end
            tau_t = 2 * pi / Nt;
            % mirror: x axis
            mirrorAxis = [1, 0];
            tmp = obj.getDefaultMeshZoneName;
            obj.changeMeshZoneName(mzname, [tmp, '1']);
            obj.cmmz([tmp, '2'], [tmp, '1'], mirrorAxis);
            
            for i = 1:2:2 * (Ns - 1)
                obj.crmz([tmp, num2str(i + 2)], [tmp, num2str(i)], tau_t);
                obj.crmz([tmp, num2str(i + 3)], [tmp, num2str(i + 1)], tau_t);
            end
            
            tmp = getlist(tmp, 1:2 * Ns);
            obj.jmzs(mzname, tmp{:});
            % change states
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function cmr(obj, mzname, Nt, Nc, type)
            tau_t = 2 * pi / Nt;
            
            switch type
                case 'x'
                    mirrorAxis = [1, 0];
                case 'y'
                    mirrorAxis = [0, 1];
                case 's'
                    mirrorAxis = [cos(tau_t / 2), sin(tau_t / 2)];
                otherwise
                    error('Wrong type.');
            end
            
            obj.changeMeshZoneName(mzname, [mzname, '11']);
            obj.cmmz([mzname, '21'], [mzname, '11'], mirrorAxis);
            
            for i = 1:(Nc - 1)
                obj.crmz([mzname, '1', num2str(i + 1)], [mzname, '1', num2str(i)], tau_t);
                obj.crmz([mzname, '2', num2str(i + 1)], [mzname, '2', num2str(i)], tau_t);
            end
            
            % change states
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function aux_cmrmzx(obj, mzname, Nt, Nc)
            % copy|mirror|rotate|x Axis
            if nargin < 4, Nc = Nt; end
            cmr(obj, mzname, Nt, Nc, 'x');
            % change states
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function aux_cmrmzs(obj, mzname, Nt, Nc)
            % copy|mirror|rotate|slot Symmetry
            if nargin < 4, Nc = Nt; end
            cmr(obj, mzname, Nt, Nc, 's');
            % change states
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function aux_crmz(obj, mzName, Nt, Nc)
            if nargin < 4, Nc = Nt; end
            tau_t = 2 * pi / Nt;
            obj.changeMeshZoneName(mzName, [mzName, '1']);
            
            for i = 1:(Nc - 1)
                obj.crmz([mzName, num2str(i + 1)], [mzName, num2str(i)], tau_t);
            end
            
            % change states
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function aux_cmyshxj(obj, mzName, Nc, xsh)
            % copy|mirror respect to Y axis|shift along X axis|join
            tmp = obj.getDefaultMeshZoneName;
            obj.changeMeshZoneName(mzName, [tmp, '1']);
            obj.cmmz([tmp, '2'], [tmp, '1'], [0, 1]);
            
            for i = 1:2:2 * (Nc - 1)
                obj.cshmz([tmp, num2str(i + 2)], [tmp, num2str(i)], [xsh, 0]);
                obj.cshmz([tmp, num2str(i + 3)], [tmp, num2str(i + 1)], [xsh, 0]);
            end
            
            tmp = getlist(tmp, 1:2 * Nc);
            obj.jmzs(mzName, tmp{:});
            % change states
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function aux_cmyshx(obj, mzName, Nc, xsh)
            % copy|mirror respect to Y axis|shift along X axis
            tmp = mzName;
            obj.changeMeshZoneName(mzName, [tmp, '11']);
            obj.cmmz([tmp, '21'], [tmp, '11'], [0, 1]);
            
            for i = 2:Nc
                obj.cshmz([tmp, '1', num2str(i)], [tmp, '1', num2str(i - 1)], [xsh, 0]);
                obj.cshmz([tmp, '2', num2str(i)], [tmp, '2', num2str(i - 1)], [xsh, 0]);
            end
            
            % change states
            obj.makeFalse_isGlobalMeshGenerated;
        end
        
        function mzptr = getRectangleEnclouser(obj, x0, y0, w, h, dx, dy)
            
            obj.ggmesh;
            eindex = find(obj.bedges);
            fi = obj.edges(obj.bedges, 1:2);
            
            for i = 1:size(fi, 1)
                
                if obj.edges(eindex(i), 4) == 0
                    fi(i, :) = fliplr(fi(i, :));
                end
                
            end
            
            p1 = obj.nodes(fi(:, 1), :);
            p2 = obj.nodes(fi(:, 2), :);
            
            Np = size(fi, 1);
            
            fi = (1:Np)';
            
            [vi, ~, index] = uniquetol([p1; p2], 1e-6, 'ByRows', true);
            
            fi = fi(index);
            fi = reshape(fi, [], 2);
            
            Nx = ceil(w / dx);
            Ny = ceil(h / dy);
            Np = 2 * Nx + 2 * (Ny - 2);
            vo = zeros(Np, 2);
            
            vo(1:Nx, 1) = linspace(x0, x0 + w, Nx);
            vo(1:Nx, 2) = y0;
            
            vo(Nx + 1:Nx + Ny - 2, 1) = x0 + w;
            tmp = linspace(y0, y0 + h, Ny);
            vo(Nx + 1:Nx + Ny - 2, 2) = tmp(2:end - 1);
            
            vo(Nx + Ny - 1:2 * Nx + Ny - 2, 1) = linspace(x0 + w, x0, Nx);
            vo(Nx + Ny - 1:2 * Nx + Ny - 2, 2) = y0 + h;
            
            vo(2 * Nx + Ny - 1:end, 1) = x0;
            tmp = linspace(y0 + h, y0, Ny);
            vo(2 * Nx + Ny - 1:end, 2) = tmp(2:end - 1);
            
            tmp = (1:Np)';
            fo = [tmp, circshift(tmp, -1)];
            
            mzptr = emdlab_m2d_mg0([fi; fo + size(vi, 1)], [vi; vo]);
            
        end

        function mzptr = getRectangleEnclouserFV(obj, x0, y0, w, h, dx, dy, f, v)
            
            obj.ggmesh;
            eindex = find(obj.bedges);
            fi = obj.edges(obj.bedges, 1:2);
            
            for i = 1:size(fi, 1)
                
                if obj.edges(eindex(i), 4) == 0
                    fi(i, :) = fliplr(fi(i, :));
                end
                
            end
            
            p1 = obj.nodes(fi(:, 1), :);
            p2 = obj.nodes(fi(:, 2), :);
            
            Np = size(fi, 1);
            
            fi = (1:Np)';
            
            [vi, ~, index] = uniquetol([p1; p2], 1e-6, 'ByRows', true);
            
            fi = fi(index);
            fi = reshape(fi, [], 2);
            
            Nx = ceil(w / dx);
            Ny = ceil(h / dy);
            Np = 2 * Nx + 2 * (Ny - 2);
            vo = zeros(Np, 2);
            
            vo(1:Nx, 1) = linspace(x0, x0 + w, Nx);
            vo(1:Nx, 2) = y0;
            
            vo(Nx + 1:Nx + Ny - 2, 1) = x0 + w;
            tmp = linspace(y0, y0 + h, Ny);
            vo(Nx + 1:Nx + Ny - 2, 2) = tmp(2:end - 1);
            
            vo(Nx + Ny - 1:2 * Nx + Ny - 2, 1) = linspace(x0 + w, x0, Nx);
            vo(Nx + Ny - 1:2 * Nx + Ny - 2, 2) = y0 + h;
            
            vo(2 * Nx + Ny - 1:end, 1) = x0;
            tmp = linspace(y0 + h, y0, Ny);
            vo(2 * Nx + Ny - 1:end, 2) = tmp(2:end - 1);
            
            tmp = (1:Np)';
            fo = [tmp, circshift(tmp, -1)];

            [f_tmp,v_tmp] = emdlab_g2d_getRectangleFVInside(x0,y0,x0+w,y0+h);
            f = [f;f_tmp+size(v,1)];
            v = [v;v_tmp];
            disFunc = @(p) ext_dpoly2d(p, f, v);
            
            mzptr = emdlab_m2d_mg10([fi; fo + size(vi, 1)], [vi; vo], disFunc);
            
        end
        
    end
    
end
