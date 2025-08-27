classdef emdlab_bcs_scalarTTL4 < handle
    
    % structures for stroing boundary conditions data
    properties (SetAccess = protected)

        dirichlet(1, 1) struct;
        oddPeriodic(1, 1) struct;
        evenPeriodic(1, 1) struct;
        neumann(1, 1) struct;
        robin(1, 1) struct;

    end
    
    % update validation
    properties (Access = protected)

        isDirichletUpdated(1, 1) logical = false;
        isOddPeriodicUpdated(1, 1) logical = false;
        isEvenPeriodicUpdated(1, 1) logical = false;
        isNeumannUpdated(1, 1) logical = false;
        isRobinUpdated(1, 1) logical = false;

    end
    
    % number of defined boundary conditions
    properties (Dependent = true, SetAccess = protected)

        Nd(1, 1) double;
        Nop(1, 1) double;
        Nep(1, 1) double;
        Nn(1, 1) double;
        Nr(1, 1) double;

    end
    
    % variables for indices and values
    properties (SetAccess = private)

        % dirichlet index
        iD(:, 1) double;

        % dirichlet value
        vD(:, 1) double;

        % odd periodic master index
        mOP(:, 1) double;

        % odd periodic slave index
        sOP(:, 1) double;

        % even periodic master index
        mEP(:, 1) double;

        % even periodic slave index
        sEP(:, 1) double;

        % index of robin edges
        iR double;

        % a value of robin edges
        aR(:, 1) double;

        % b value of robin edges
        bR(:, 1) double;

        % index of neumann edges
        iN double;

        % b value of robin edges
        bN(:, 1) double;

    end
    
    % number of indices values
    properties (Dependent = true, SetAccess = protected)

        % number of diriclet nodes
        Ndbcs(1, 1) double;

        % number of odd periodic nodes
        Nopbcs(1, 1) double;

        % number of even periodic nodes
        Nepbcs(1, 1) double;

        % number of neumann edges or facets
        Nnbcs(1, 1) double;
        
        % number of robin edges or facets
        Nrbcs(1, 1) double;
        
    end
    
    methods
        %% constructor
        function obj = emdlab_bcs_scalarTTL4()
        end
        
        %% Editing Functions
        function y = getBCNames(obj)
            y = [fieldnames(obj.dirichlet); ...
                fieldnames(obj.oddPeriodic); ...
                fieldnames(obj.evenPeriodic); ...
                fieldnames(obj.neumann); ...
                fieldnames(obj.robin)];
        end
        
        function y = getDefaultBCName(obj)
            index = 0;
            bcNames = obj.getBCNames;
            
            while true
                index = index + 1;
                y = ['BC', num2str(index)];
                if ~ismember(y, bcNames), break; end
            end
            
        end
        
        function bcName = checkBCExictence(obj, bcName)
            bcName = rmspaces(bcName);
            
            if ~ismember(bcName, obj.getBCNames)
                throw(MException('', ['Boundary condition with the name <', bcName, '> does not exist.']));
            end
            
        end
        
        function bcName = checkBCNonExictence(obj, bcName)
            bcName = rmspaces(bcName);
            
            if ismember(bcName, obj.getBCNames)
                throw(MException('', ['Boundary condition with the name <', bcName, '> already exist.']));
            end
            
        end
        
        function clearAllBCs(obj)
            % removing structurs
            obj.dirichlet = struct();
            obj.oddPeriodic = struct();
            obj.evenPeriodic = struct();
            obj.neumann = struct();
            obj.robin = struct();
            % changing updates
            obj.isDirichletUpdated = false;
            obj.isNeumannUpdated = false;
            obj.isRobinUpdated = false;
            obj.isOddPeriodicUpdated = false;
            obj.isEvenPeriodicUpdated = false;
        end
        
        %% getting number of defined bouandary conditions
        function y = get.Nd(obj)
            y = numel(fieldnames(obj.dirichlet));
        end
        
        function y = get.Nop(obj)
            y = numel(fieldnames(obj.oddPeriodic));
        end
        
        function y = get.Nep(obj)
            y = numel(fieldnames(obj.evenPeriodic));
        end
        
        function y = get.Nn(obj)
            y = numel(fieldnames(obj.neumann));
        end
        
        function y = get.Nr(obj)
            y = numel(fieldnames(obj.robin));
        end
        
        %% definitions
        function setDirichlet(obj, index, value, bcName)
            
            if nargin < 4
                bcName = obj.getDefaultBCName;
            end
            
            bcName = obj.checkBCNonExictence(bcName);
            [m, n] = size(value);
            
            if (m == 1 && n == 1) || (m == size(index, 1) && n == 1)
                obj.dirichlet.(bcName).index = index;
                obj.dirichlet.(bcName).value = value;
            else
                throw(MException('', 'Value is in wrong format.'));
            end
            
            % change states
            obj.isDirichletUpdated = false;
        end
        
        function setRobin(obj, index, alpha, beta, bcName)
            
            if nargin < 5
                bcName = obj.getDefaultBCName;
            end
            
            bcName = obj.checkBCNonExictence(bcName);
            [mAlpha, nAlpha] = size(alpha);
            [mBeta, nBeta] = size(alpha);
            
            if (mAlpha == 1 && nAlpha == 1 && mBeta == 1 && nBeta == 1)
                obj.robin.(bcName).index = index;
                obj.robin.(bcName).alpha = alpha;
                obj.robin.(bcName).beta = beta;
            else
                throw(MException('', 'Value is in wrong format.'));
            end
            
            % change states
            obj.isRobinUpdated = false;
        end
        
        function setNeumann(obj, index, beta, bcName)
            
            if nargin < 4
                bcName = obj.getDefaultBCName;
            end
            
            bcName = obj.checkBCNonExictence(bcName);
            [mBeta, nBeta] = size(beta);
            
            if (mBeta == 1 && nBeta == 1)
                obj.neumann.(bcName).index = index;
                obj.neumann.(bcName).beta = beta;
            else
                throw(MException('', 'Value is in wrong format.'));
            end
            
            % change states
            obj.isNeumannUpdated = false;
        end
        
        function setOddPeriodic(obj, mIndex, sIndex, bcName)
            
            if nargin < 4
                bcName = obj.getDefaultBCName;
            end
            
            bcName = obj.checkBCNonExictence(bcName);
            [mm, nm] = size(mIndex);
            [ms, ns] = size(sIndex);
            
            if (mm == ms && nm == 1 && ns == 1)
                obj.oddPeriodic.(bcName).mIndex = mIndex;
                obj.oddPeriodic.(bcName).sIndex = sIndex;
            else
                throw(MException('', 'Invalid mIndex and sIndex format.'));
            end
            
            % change states
            obj.isOddPeriodicUpdated = false;
        end
        
        function setEvenPeriodic(obj, mIndex, sIndex, bcName)
            
            if nargin < 4
                bcName = obj.getDefaultBCName;
            end
            
            bcName = obj.checkBCNonExictence(bcName);
            [mm, nm] = size(mIndex);
            [ms, ns] = size(sIndex);
            
            if (mm == ms && nm == 1 && ns == 1)
                obj.evenPeriodic.(bcName).mIndex = mIndex;
                obj.evenPeriodic.(bcName).sIndex = sIndex;
            else
                throw(MException('', 'Invalid mIndex and sIndex format.'));
            end
            
            % change states
            obj.isEvenPeriodicUpdated = false;
        end
        
        %% getters
        function y = get.Ndbcs(obj)
            y = length(obj.iD);
        end
        
        function y = get.Nopbcs(obj)
            y = length(obj.mOP);
        end
        
        function y = get.Nepbcs(obj)
            y = length(obj.mEP);
        end
        
        function y = get.Nrbcs(obj)
            y = length(obj.iR);
        end
        
        function y = get.Nnbcs(obj)
            y = length(obj.iN);
        end
        
        %% updates
        function updateDirichlet(obj)
            dBCs = fieldnames(obj.dirichlet);
            NdBCs = numel(dBCs);
            tmp = zeros(1, NdBCs + 1);
            
            for i = 1:NdBCs
                tmp(i + 1) = length(obj.dirichlet.(dBCs{i}).index);
            end
            
            N = sum(tmp);
            
            for i = 2:length(tmp)
                tmp(i) = tmp(i) + tmp(i - 1);
            end
            
            obj.iD = zeros(N, 1);
            obj.vD = zeros(N, 1);
            
            for i = 1:NdBCs
                obj.iD(tmp(i) + 1:tmp(i + 1)) = obj.dirichlet.(dBCs{i}).index;
                obj.vD(tmp(i) + 1:tmp(i + 1)) = obj.dirichlet.(dBCs{i}).value;
            end
            
        end
        
        function updateNeumann(obj)
            nBCs = fieldnames(obj.neumann);
            NnBCs = numel(nBCs);
            tmp = zeros(1, NnBCs + 1);
            
            for i = 1:NnBCs
                tmp(i + 1) = size(obj.neumann.(nBCs{i}).index, 1);
            end
            
            N = sum(tmp);
            
            for i = 2:length(tmp)
                tmp(i) = tmp(i) + tmp(i - 1);
            end
            
            switch obj.eType
                case 'TL3'
                    Dim = 2;
                case 'TTL4'
                    Dim = 3;
                otherwise
                    throw(MException('', 'Wrong element type.'));
            end
            
            obj.iN = zeros(N, Dim);
            obj.bN = zeros(N, 1);
            
            for i = 1:NnBCs
                bcptr = obj.neumann.(nBCs{i});
                obj.iN(tmp(i) + 1:tmp(i + 1), :) = bcptr.index;
                obj.bN(tmp(i) + 1:tmp(i + 1)) = bcptr.beta;
            end
            
        end
        
        function updateRobin(obj)
            rBCs = fieldnames(obj.robin);
            NrBCs = numel(rBCs);
            tmp = zeros(1, NrBCs + 1);
            
            for i = 1:NrBCs
                tmp(i + 1) = size(obj.robin.(rBCs{i}).index, 1);
            end
            
            N = sum(tmp);
            
            for i = 2:length(tmp)
                tmp(i) = tmp(i) + tmp(i - 1);
            end
            
            switch obj.eType
                case 'TL3'
                    Dim = 2;
                case 'TTL4'
                    Dim = 3;
                otherwise
                    throw(MException('', 'Wrong element type.'));
            end
            
            obj.iR = zeros(N, Dim);
            obj.aR = zeros(N, 1);
            obj.bR = zeros(N, 1);
            
            for i = 1:NrBCs
                bcptr = obj.robin.(rBCs{i});
                obj.iR(tmp(i) + 1:tmp(i + 1), :) = bcptr.index;
                obj.aR(tmp(i) + 1:tmp(i + 1)) = bcptr.alpha;
                obj.bR(tmp(i) + 1:tmp(i + 1)) = bcptr.beta;
            end
            
        end
        
        function updateOddPeriodic(obj)
            opBCs = fieldnames(obj.oddPeriodic);
            NopBCs = numel(opBCs);
            tmp = zeros(1, NopBCs + 1);
            
            for i = 1:NopBCs
                tmp(i + 1) = length(obj.oddPeriodic.(opBCs{i}).mIndex);
            end
            
            N = sum(tmp);
            
            for i = 2:length(tmp)
                tmp(i) = tmp(i) + tmp(i - 1);
            end
            
            obj.mOP = zeros(N, 1);
            obj.sOP = zeros(N, 1);
            
            for i = 1:NopBCs
                obj.mOP(tmp(i) + 1:tmp(i + 1)) = obj.oddPeriodic.(opBCs{i}).mIndex;
                obj.sOP(tmp(i) + 1:tmp(i + 1)) = obj.oddPeriodic.(opBCs{i}).sIndex;
            end
            
        end
        
        function updateEvenPeriodic(obj)
            epBCs = fieldnames(obj.evenPeriodic);
            NepBCs = numel(epBCs);
            tmp = zeros(1, NepBCs + 1);
            
            for i = 1:NepBCs
                tmp(i + 1) = length(obj.evenPeriodic.(epBCs{i}).mIndex);
            end
            
            N = sum(tmp);
            
            for i = 2:length(tmp)
                tmp(i) = tmp(i) + tmp(i - 1);
            end
            
            obj.mEP = zeros(N, 1);
            obj.sEP = zeros(N, 1);
            
            for i = 1:NepBCs
                obj.mEP(tmp(i) + 1:tmp(i + 1)) = obj.evenPeriodic.(epBCs{i}).mIndex;
                obj.sEP(tmp(i) + 1:tmp(i + 1)) = obj.evenPeriodic.(epBCs{i}).sIndex;
            end
            
        end
        
        function updateAll(obj)
            obj.updateDirichlet;
            obj.updateNeumann;
            obj.updateRobin;
            obj.updateOddPeriodic;
            obj.updateEvenPeriodic;
        end
        
        %% imposing of boundary conditions
        function [K, F] = imposeDirichlet(obj, K, F)
            
            if obj.Nd
                F(obj.iD) = obj.vD;
                K(obj.iD, :) = sparse(1:obj.Ndbcs, obj.iD, ones(1, obj.Ndbcs), obj.Ndbcs, length(F));
            end
            
        end
        
        function [K, F] = imposeZeroDirichlet(obj, K, F)
            
            if obj.Nd
                F(obj.iD) = 0;
                K(obj.iD, :) = sparse(1:obj.Ndbcs, obj.iD, ones(1, obj.Ndbcs), obj.Ndbcs, length(F));
            end
            
        end
        
    end
    
end
