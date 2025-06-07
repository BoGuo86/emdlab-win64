% developer: https://ComProgExpert.com
% loop class for 2D geometries
% loop is a sequence of edges, which is alway counterclockwise

classdef emdlab_g2d_loop < handle

    properties

        % a cell to store edges
        edges (1,:) cell;

        % direction of each edge
        directions (1,:) logical;

        % the object tag
        tag (1,:) char = '';

    end

    methods

        % constructor
        function  obj = emdlab_g2d_loop()
        end

        % add a new edge
        function addEdge(obj, newEdge, direction)

            obj.edges{end+1} = newEdge;
            obj.directions(end+1) = direction;

        end

        % clear all existing edges
        function clearAllEdges(obj)

            obj.edges = {};
            obj.directions = [];

        end

        % get mesh nodes of a loop
        function [nodes, cl] = getMeshNodes(obj)

            % number of edges
            Ne = numel(obj.edges);

            % store edge nodes and edge connectivity list in cells
            nodes = cell(Ne,1);

            % get edges mesh nodes
            for i = 1:Ne
                nodes{i} = obj.edges{i}.getMeshNodes;
                if ~obj.directions(i)
                    nodes{i} = flipud(nodes{i});
                end
                nodes{i}(end,:) = [];
            end

            % unify nodes
            nodes = cell2mat(nodes);

            % make the connectivity list
            cl = (1:size(nodes,1))';
            cl = [cl, circshift(cl,-1)];

        end

        function [nodes, cl] = getMeshNodesMinimal(obj)

            % number of edges
            Ne = numel(obj.edges);

            % store edge nodes and edge connectivity list in cells
            nodes = cell(Ne,1);

            % get edges mesh nodes
            for i = 1:Ne
                nodes{i} = obj.edges{i}.getMeshNodesMinimal;
                if ~obj.directions(i)
                    nodes{i} = flipud(nodes{i});
                end
                nodes{i}(end,:) = [];
            end

            % unify nodes
            nodes = cell2mat(nodes);

            % make the connectivity list
            cl = (1:size(nodes,1))';
            cl = [cl, circshift(cl,-1)];

        end

    end

end