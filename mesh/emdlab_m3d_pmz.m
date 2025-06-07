classdef emdlab_m3d_pmz < handle & matlab.mixin.Copyable & emdlab_g2d_constants
    
    properties (SetAccess = protected)
        % mesh nodes
        nodes(3, :) double;
        % mesh connectivity list
        cl(6, :) double;
        % mesh elements
        elements(5, :) double;
        % tfacets
        tfacets(3, :) double;
        % qfacets
        qfacets(4, :) double;
        % edges
        edges(2,:) double;
        % Named Selections
        nodeNamedSelections(1, 1) struct;
    end
    
    properties (Access = private)
        % element area
        ea(:, 1) double;
        % mesh zone area
        area(1, 1) double;
        % Q
        Q(1, :) double;
        % Weight matrix
        Wm
        % states
        isDataSetted(1, 1) logical = false;
        is_ea_Evaluated(1, 1) logical = false;
        is_area_Evaluated(1, 1) logical = false;
        is_Wm_Evaluated(1, 1) logical = false;
        is_Q_Evaluated(1, 1) logical = false;
    end
    
    properties (Dependent = true)
        % number of mesh zone nodes
        Nn(1, 1) double;
        % number of mesh zone elements
        Ne(1, 1) double;
    end
    
    properties
        % zone index
        zi(1, 1) double;
        % local to global node index
        l2g(:, 1) double;
        % material of zone
        material char = 'air';
        % mesh zone color
        color = 'c';
        % mesh zone properties: differs in differents solvers
        props(1, 1) struct;
    end
    
    methods
        
        function obj = emdlab_m3d_pmz(cl, nodes)
            obj.cl = cl;
            obj.nodes = nodes;
        end
        
        function y = get.Nn(obj)
            y = size(obj.nodes, 2);
        end
        
        function y = get.Ne(obj)
            y = size(obj.cl, 2);
        end
        
        function setdata(obj)
            
            
            e1 = obj.cl([1, 2], :);
            e2 = obj.cl([2, 3], :);
            e3 = obj.cl([3, 1], :);
            e4 = obj.cl([4, 5], :);
            e5 = obj.cl([5, 6], :);
            e6 = obj.cl([6, 4], :);
            e7 = obj.cl([1, 4], :);
            e8 = obj.cl([2, 5], :);
            e9 = obj.cl([3, 6], :);
            
            
            % sorting for lower index
            [e1, s1] = sort(e1);
            [e2, s2] = sort(e2);
            [e3, s3] = sort(e3);
            [e4, s4] = sort(e4);
            [e5, s5] = sort(e5);
            [e6, s6] = sort(e6);
            [e7, s7] = sort(e7);
            [e8, s8] = sort(e8);
            [e9, s9] = sort(e9);
            
            % specefying changed edge index
            s1 = s1(1,:) == 2;
            s2 = s2(1,:) == 2;
            s3 = s3(1,:) == 2;
            s4 = s4(1,:) == 2;
            s5 = s5(1,:) == 2;
            s6 = s6(1,:) == 2;
            s7 = s7(1,:) == 2;
            s8 = s8(1,:) == 2;
            s9 = s9(1,:) == 2;
            
            % unification of edges
            [tmp, ~, ic] = unique([e1, e2, e3, e4, e5, e6, e7, e8, e9]', 'rows');
            
            obj.edges = tmp';
            % getting number of elements
            ne = obj.Ne;
            % getting index of edge corresponding to each elements
            e1 = ic(1:ne);
            e2 = ic(1 + ne:2 * ne);
            e3 = ic(1 + 2 * ne:3 * ne);
            e4 = ic(1 + 3 * ne:4 * ne);
            e5 = ic(1 + 4 * ne:5 * ne);
            e6 = ic(1 + 5 * ne:6 * ne);
            e7 = ic(1 + 6 * ne:7 * ne);
            e8 = ic(1 + 7 * ne:8 * ne);
            e9 = ic(1 + 8 * ne:9 * ne);
            
            e1(s1) = -e1(s1);
            e2(s2) = -e2(s2);
            e3(s3) = -e3(s3);
            e4(s4) = -e4(s4);
            e5(s5) = -e5(s5);
            e6(s6) = -e6(s6);
            e7(s7) = -e7(s7);
            e8(s8) = -e8(s8);
            e9(s9) = -e9(s9);
            
            obj.tfacets = [e1,e2,e3;e4,e5,e6]';
            obj.qfacets = [e1,e8,e4,e7;e2,e9,e5,e8;e3,e7,e6,e9]';
            
        end
    end
    
end
