classdef emdlab_solvers_IHLMSTL6 < handle
    properties
        % solver mesh
        m
        % boundary conditions
        dbcs
        opbcs
        epbcs
        % units class
        units = UNITSC;
        % physical constants
        pcts = PCONSTANTS;
        % magnetic vector potential: goal quantity
        A
        % magnetic field density or magnetic induction
        Be3in
        Be3out
        Be6
        % nodal B
        Bn
        % magnetic field intensity
        H
        % elements data
        edata
    end
    properties (SetAccess = private)
        % depth of problem
        depth = 1000;
    end
    properties (Access = private)
        isBnEvaluated = false;
        isElementDataAssigned = false;
        is_Be3in_Evaluated = false;
        is_Be3out_Evaluated = false;
        is_Be6_Evaluated = false;
    end
    methods
        %% Initialization
        function obj = emdlab_solvers_IHLMSTL6(m)
            obj.m = m;
            % default values of bcs
            obj.dbcs = zeros([],2);
            obj.opbcs = zeros([],2);
            obj.epbcs = zeros([],2);
            % default properties of mzs
            mznames = fieldnames(obj.m.mzs);
            for i = 1:numel(mznames)
                obj.setdp(mznames{i});
            end
        end
        function set.m(obj,m)
            if ~isa(m,'TMDBC')
                error('Input must be a mesh class of type TMDBC.');
            end
            obj.m = m;
        end
        function setDepth(obj,value)
            obj.depth = value;
        end
        %% Solver Properties for Mesh Zones
        function setdp(obj,mzname)
            obj.m.mzs.(mzname).props.isExcited = false;
            obj.m.mzs.(mzname).props.isMagnetized = false;
            obj.makeFalse_isElementDataAssigned;
        end
        function addmz(obj,mzname,mzvalue)
            mzname = rmspaces(mzname);
            obj.m.addmz(mzname,mzvalue);
            obj.setdp(mzname);
            obj.makeFalse_isElementDataAssigned;
        end
        function removemz(obj,mzname)
            obj.m.removemz(mzname);
            obj.makeFalse_isElementDataAssigned;
        end
        %% Setting Excitations
        function setExcitation(obj,mzname,value,varargin)
            % varargin = Excitation type: Current or Current Density
            mzname = rmspaces(mzname);
            mznames = fieldnames(obj.m.mzs);
            if ~ismember(mzname,mznames)
                error('Specefied zone does not exist.');
            end
            obj.m.mzs.(mzname).props.isExcited = true;
            obj.m.mzs.(mzname).props.excitation = msExcitation(value,varargin{:});
            % change states
            obj.makeFalse_isElementDataAssigned;
        end
        function setMagnetization(obj,mzname,value)
            mzname = rmspaces(mzname);
            mznames = fieldnames(obj.m.mzs);
            if ~ismember(mzname,mznames)
                error('Specefied zone does not exist.');
            end
            if ~(isa(value,'function_handle')||(isa(value,'msMagnetization')))
                error('Magnetization must be a function handle.');
            end
            obj.m.mzs.(mzname).props.isMagnetized = true;
            obj.m.mzs.(mzname).props.magnetization = value;
            % change states
            obj.makeFalse_isElementDataAssigned;
        end
        %% Boundary Conditions
        function setdbc(obj,k,value)
            Nk = length(k);
            Nvalue = length(value);
            if Nvalue == 1
                obj.dbcs = [obj.dbcs;k,value*ones(Nk,1)];
            elseif Nk == Nvalue
                obj.dbcs = [obj.dbcs;k,value];
            else
                error('Value mest be scalar or a vector with the same length of k.');
            end
        end
        function setopbc(obj,km,ks)
            if length(km) == length(ks)
                obj.opbcs = [obj.opbcs;km,ks];
            else
                error('The length of vectors km and ks must be the same.');
            end
        end
        function setepbc(obj,km,ks)
            if length(km) == length(ks)
                obj.epbcs = [obj.epbcs;km,ks];
            else
                error('The length of vectors km and ks must be the same.');
            end
        end
        function cleardbcs(obj)
            obj.dbcs = zeros([],2);
        end
        function clearopbcs(obj)
            obj.opbcs = zeros([],2);
        end
        function clearepbcs(obj)
            obj.epbcs = zeros([],2);
        end
        function clearallbcs(obj)
            obj.cleardbcs;
            obj.clearopbcs;
            obj.clearepbcs;
        end
        %% Solver
        function obj = solve(obj)
            obj.assignEdata;
            % Construction of [K] and [F]
            tic
            disp('*******************************************************')
            % Assembeling [F]
%             fi = repmat(obj.edata.InternalCurrentDensity,6,1);
%             Mx = repmat(obj.edata.Magnetization(1,:),3,1)*sparse(1:obj.m.Ne,1:obj.m.Ne,obj.m.gea);
%             My = repmat(obj.edata.Magnetization(2,:),3,1)*sparse(1:obj.m.Ne,1:obj.m.Ne,obj.m.gea);
%             F = (fi.*obj.m.Fe*obj.units.K_currentDensity+...
%                 (obj.m.gphiy.*Mx-obj.m.gphix.*My)/obj.units.K_length);
% Assembeling [F]
            F = sparse(obj.m.cl',ones(6*obj.m.Ne,1),repmat(obj.edata.InternalCurrentDensity,6,1).*obj.m.Fe);
            % applying scales on load vector
            F = F*obj.units.K_currentDensity*obj.units.K_length^2/obj.units.K_magneticVectorPotential;
%             F = sparse(obj.m.cl',ones(3*obj.m.Ne,1),F);
            % Assembeling [K]
            [Iindex,Jindex] = getij(6,1);
            K = sparse(obj.m.cl(:,Iindex)',...
                obj.m.cl(:,Jindex)',...
                repmat(obj.edata.MagneticReluctivity,36,1).* ...
                obj.m.Ke(getkindex(6),:));
            disp('Construction of [K] and [F] compeleted.')
            toc
            % imposing boundary conditions on [K] and [F]
            % dbcs
            if ~isempty(obj.dbcs)
                Ndbcs = length(obj.dbcs(:,1));
                F(obj.dbcs(:,1)) = obj.dbcs(:,2);
                K(obj.dbcs(:,1),:) = sparse(1:Ndbcs,...
                    obj.dbcs(:,1),ones(1,Ndbcs),Ndbcs,obj.m.Nn);
            end
            % opbcs
            if ~isempty(obj.opbcs)
                Nopbcs = size(obj.opbcs,1);
                F(obj.opbcs(:,1)) = F(obj.opbcs(:,1)) - F(obj.opbcs(:,2));
                F(obj.opbcs(:,2)) = 0;
                K(obj.opbcs(:,1),:) = K(obj.opbcs(:,1),:) - K(obj.opbcs(:,2),:);
                K(obj.opbcs(:,2),:) = sparse([1:Nopbcs,1:Nopbcs],...
                    obj.opbcs(:),ones(1,2*Nopbcs),Nopbcs,obj.m.Nn);
            end
            % epbcs
            if ~isempty(obj.epbcs)
                Nepbcs = size(obj.epbcs,1);
                F(obj.epbcs(:,1)) = F(obj.epbcs(:,1)) + F(obj.epbcs(:,2));
                F(obj.epbcs(:,2)) = 0;
                K(obj.epbcs(:,1),:) = K(obj.epbcs(:,1),:) + K(obj.epbcs(:,2),:);
                K(obj.epbcs(:,2),:) = sparse([1:Nepbcs,1:Nepbcs],...
                    obj.epbcs(:),[ones(1,Nepbcs),-ones(1,Nepbcs)],Nepbcs,obj.m.Nn);
            end
            disp('All boundary condition imposed')
            toc
            % solving [K][U] = [F]
            tic
            disp('*******************************************************')
            % solving equation KU = F
            if ~any(F)
                obj.A = full(F);
                return
            end
            obj.A = full(K\F);
%             obj.evalB;
            disp('initial geuss evaluated.')
            toc
        end
        function assignEdata(obj)
            if obj.isElementDataAssigned
                return
            end
            % preparing mesh data
            obj.m.ggmesh;
            obj.m.gd2elements;
            obj.m.evalKeFe_TL6;
            tic
            % assigning material and force data to each triangle
            % initialization
            obj.edata.MagneticReluctivity = zeros(1,obj.m.Ne);
            obj.edata.ElectricConductivity = zeros(1,obj.m.Ne);
            obj.edata.InternalCurrentDensity = zeros(1,obj.m.Ne);
            obj.edata.Magnetization = zeros(2,obj.m.Ne);
            mznames = fieldnames(obj.m.mzs);
            for i = 1:obj.m.Nmzs
                mzname = mznames{i};
                if ~obj.m.mts.(obj.m.mzs.(mzname).material).MagneticPermeability.isIsotropic
                    error('Some materials are NonIsotropic.');
                elseif ~obj.m.mts.(obj.m.mzs.(mzname).material).MagneticPermeability.isLinear
                    error('Some materials are Non-Linear.');
                else
                    % assigning Magnetic Permeability
                    obj.edata.MagneticReluctivity(obj.m.ezi(:,obj.m.mzs.(mzname).zi)) =...
                        1/obj.m.mts.(obj.m.mzs.(mzname).material).MagneticPermeability.value;
                end
                % assigning Electric Conductivity
                obj.edata.ElectricConductivity(obj.m.ezi(:,obj.m.mzs.(mzname).zi)) =...
                    obj.m.mts.(obj.m.mzs.(mzname).material).ElectricConductivity.value;
                % assigning Internal Current Density
                if obj.m.mzs.(mzname).props.isExcited
                    switch obj.m.mzs.(mzname).props.excitation.type
                        case 'currentDensity'
                            obj.edata.InternalCurrentDensity(obj.m.ezi(:,obj.m.mzs.(mzname).zi)) = ...
                                obj.m.mzs.(mzname).props.excitation.value;
                        case 'current'
                            obj.edata.InternalCurrentDensity(obj.m.ezi(:,obj.m.mzs.(mzname).zi)) = ...
                                obj.m.mzs.(mzname).props.excitation.value/obj.m.mzs.(mzname).area;
                    end
                end
                % assigning Magnetization
                if obj.m.mzs.(mzname).props.isMagnetized
                    index = obj.m.cl(obj.m.ezi(:,obj.m.mzs.(mzname).zi),:);
                    centers = (obj.m.nodes(index(:,1),:)+obj.m.nodes(index(:,2),:)+obj.m.nodes(index(:,3),:))/3;
                    if isa(obj.m.mzs.(mzname).props.magnetization,'msMagnetization')
                        M = obj.m.mzs.(mzname).props.magnetization.getM(centers);
                    else
                        M = feval(obj.m.mzs.(mzname).magnetization,centers);
                    end
                    obj.edata.Magnetization(1,obj.m.ezi(:,obj.m.mzs.(mzname).zi)) = M(:,1)';
                    obj.edata.Magnetization(2,obj.m.ezi(:,obj.m.mzs.(mzname).zi)) = M(:,2)';
                    
                end
            end
            disp('Initialization of material and force data compeleted.')
            toc
            % change states
            obj.isElementDataAssigned = true;
        end
        %% PostProccessing
        % evaluations
        function evalBe3out(obj)
            if obj.is_Be3out_Evaluated
                return;
            end  
            [obj.Be3out.x,obj.Be3out.y] = IHLMSTL6_evalBe3out(obj.m.cl,obj.A,....
                obj.m.Aux.JITGG00x,obj.m.Aux.JITGG00y,...
                obj.m.Aux.JITGG10x,obj.m.Aux.JITGG10y,...
                obj.m.Aux.JITGG01x,obj.m.Aux.JITGG01y);
            obj.Be3out.x = obj.Be3out.x *(obj.units.K_magneticVectorPotential...
                /obj.units.K_length);
            obj.Be3out.y = obj.Be3out.y *(obj.units.K_magneticVectorPotential...
                /obj.units.K_length);
            % change states
            obj.is_Be3out_Evaluated = true;
        end
        function evalBe3in(obj)
            if obj.is_Be3in_Evaluated
                return;
            end  
            [obj.Be3in.x,obj.Be3in.y] = IHLMSTL6_evalBe3in(obj.m.cl,obj.A,....
                obj.m.Aux.JITGG120x,obj.m.Aux.JITGG120y,...
                obj.m.Aux.JITGG1212x,obj.m.Aux.JITGG1212y,...
                obj.m.Aux.JITGG012x,obj.m.Aux.JITGG012y);
            obj.Be3in.x = obj.Be3in.x *(obj.units.K_magneticVectorPotential...
                /obj.units.K_length);
            obj.Be3in.y = obj.Be3in.y *(obj.units.K_magneticVectorPotential...
                /obj.units.K_length);
            % change states
            obj.is_Be3in_Evaluated = true;
        end
		function evalBe6(obj)
            if obj.is_Be6_Evaluated
                return;
            end
            [obj.Be6.x,obj.Be6.y] = IHLMSTL6_evalBe6(obj.m.cl,obj.A,....
                obj.m.Aux.JITGG00x,obj.m.Aux.JITGG00y,...
                obj.m.Aux.JITGG10x,obj.m.Aux.JITGG10y,...
                obj.m.Aux.JITGG01x,obj.m.Aux.JITGG01y,...
                obj.m.Aux.JITGG120x,obj.m.Aux.JITGG120y,...
                obj.m.Aux.JITGG1212x,obj.m.Aux.JITGG1212y,...
                obj.m.Aux.JITGG012x,obj.m.Aux.JITGG012y);
            obj.Be6.x = obj.Be6.x *(obj.units.K_magneticVectorPotential...
                /obj.units.K_length);
            obj.Be6.y = obj.Be6.y *(obj.units.K_magneticVectorPotential...
                /obj.units.K_length);
            % change states
            obj.is_Be6_Evaluated = true;
        end
        function y = getBeMean(obj)
            obj.evalBe3out;
            y = [mean(obj.Be3out.x,2),mean(obj.Be3out.y,2)];
        end
        function evalBn(obj)
            if obj.isBnEvaluated
                return
            end
            Wm = sparse(obj.m.cl(:),repmat((1:obj.m.Ne)',3,1),ones(3*obj.m.Ne,1));
            obj.Bn = zeros(obj.m.Nn,2);
            obj.Bn(:,1) = (Wm*(obj.B(:,1).*obj.m.gea'))./(Wm*obj.m.gea');
            obj.Bn(:,2) = (Wm*(obj.B(:,2).*obj.m.gea'))./(Wm*obj.m.gea');
            obj.isBnEvaluated = true;
        end
        function y = evalTotalEnergy(obj)
            y = obj.units.K_length^2*0.5*...
                ((obj.m.gea.*obj.edata.MagneticReluctivity)*...
                sum(obj.B.*obj.B,2));
        end
        function y = evalLF(obj,mzname)
            mzname = rmspaces(mzname);
            y = obj.m.gea(obj.m.ezi(:,obj.m.mzs.(mzname).zi))*sum(obj.A(obj.m.cl(obj.m.ezi(:,obj.m.mzs.(mzname).zi),4:6)),2);
            y = y*obj.getDepth/3/obj.m.mzs.(mzname).area;
        end
        function plotBmag(obj,varargin)
            obj.evalBe3out;
            % specefying zones
            if ~numel(varargin)
                ti = true(obj.m.Ne,1);
            else
                ti = false(obj.m.Ne,1);
                for i = 1:numel(varargin)
                    ti = bitor(ti,obj.m.ezi(:,obj.m.mzs.(rmspaces(varargin{i})).zi));
                end
            end
            axis off equal;
            ampB = sqrt(obj.Be3out.x.^2+obj.Be3out.y.^2)'; 
            ampB = ampB(:);
            xn = obj.m.nodes(:,1);
            yn = obj.m.nodes(:,2);
            xn = xn(obj.m.cl(:,1:3)');
            yn = yn(obj.m.cl(:,1:3)');
            tmp = (1:3*obj.m.Ne)';
            tmp = reshape(tmp,3,[]);
            tmp = tmp';
            % plot through patch
            patch('Faces',tmp,'Vertices',[xn(:),yn(:)],...
                'FaceVertexCData',ampB,'FaceColor','interp',...
                'EdgeColor','none')
            view([0,0,1]);
            c = colorbar;
            c.Color = 'w';
            c.Ticks = linspace(min(ampB),max(ampB),10);
            c.Limits = [min(ampB),max(ampB)];
            c.FontSize = 20;
            colormap jet;
            title('Magnetic Flux Density Amplitude [Tesla]',...
                'color','w','fontsize',15);
            hold on
            index = obj.m.edges(:,3) - obj.m.edges(:,4);
            patch('faces',obj.m.edges(logical(abs(index)),1:2),'vertices',obj.m.nodes,...
                'EdgeColor','w');
            set(gcf,'Color',[0.2353 0.2353 0.2353])
            set(gcf,'Renderer','opengl');
            set(gca,'Clipping','off');
            zoom on;
        end
        function plotBmagSmooth(obj,varargin)
            % specefying zones
            if ~numel(varargin)
                ti = true(obj.m.Ne,1);
            else
                ti = false(obj.m.Ne,1);
                for i = 1:numel(varargin)
                    ti = bitor(ti,obj.m.ezi(:,obj.m.mzs.(rmspaces(varargin{i})).zi));
                end
            end
            f = figure('Name','Bmag','NumberTitle','off',...
                'NextPlot','new');
            a = axes(f);
            axis off equal;
            obj.evalBn;
            ampB = sqrt(sum(obj.Bn.^2,2));
            % plot through patch
            patch('Faces',obj.m.cl(ti,1:3),'Vertices',obj.m.nodes,...
                'FaceVertexCData',ampB,'FaceColor','interp',...
                'EdgeColor','none','Parent',a)
            view([0,0,1]);
            make_colorbar(min(ampB),max(ampB));
            colormap jet;
            set(f,'Color','k')
            set(gca,'Clipping','off')
            title('Magnetic Flux Density Amplitude [Tesla]',...
                'color','c','fontsize',15);
            set(f,'Renderer','opengl');
            hold on
            index = obj.m.edges(:,3) - obj.m.edges(:,4);
            patch('faces',obj.m.edges(logical(abs(index)),1:2),'vertices',obj.m.nodes,...
                'EdgeColor','w','Parent',a);
            zoom on;
        end
        function plotBvec(obj,varargin)
            if ~numel(varargin)
                ti = true(obj.m.Ne,1);
            else
                ti = false(obj.m.Ne,1);
                for i = 1:numel(varargin)
                    ti = bitor(ti,obj.m.ezi(:,obj.m.mzs.(rmspaces(varargin{i})).zi));
                end
            end
            tr = triangulation(obj.m.cl(:,1:3),obj.m.nodes);
            c = tr.incenter;
            f = figure('Name','Bvec','NumberTitle','off',...
                'NextPlot','new','WindowStyle','modal');
            axis off equal
            a = gca;
            axis off equal; hold on;
            set(f,'Color',[0.2353 0.2353 0.2353])
            quiver(c(ti,1),c(ti,2),obj.B(ti,1),obj.B(ti,2),...
                'color','w','parent',a);
            title('Magnetic Flux Density Vector','color','c','fontsize',15);
            hold on
            index = obj.m.edges(:,3) - obj.m.edges(:,4);
            patch('faces',obj.m.edges(logical(abs(index)),1:2),'vertices',obj.m.nodes,...
                'EdgeColor','w','Parent',a);
            set(f,'Renderer','opengl');
            zoom on;
        end
        function plotMvec(obj,scale)
            if nargin<2
                scale = 1;
            end
            tr = triangulation(obj.m.cl(:,1:3),obj.m.nodes);
            c = tr.incenter;
            close all;
            axis off equal; hold on;
            set(gcf,'Color','k')
            quiver(c(:,1),c(:,2),obj.edata.Magnetization(1,:)',obj.edata.Magnetization(2,:)',...
                scale,'color','w');
            title('Magnetization Vector','color','c','fontsize',15);
            set(gcf,'Renderer','opengl');
            zoom on; hold on
        end 
        function plotAmag1(obj)
            f = figure('Name','Amag','NumberTitle','off',...
                'NextPlot','new','WindowStyle','modal','Renderer','opengl',...
                'Color','k');
            a = axes(f);
            axis off equal
            patch('faces',obj.m.cl(:,1:3),'vertices',obj.m.nodes,...
                'FaceVertexCData',obj.A,'FaceColor','interp',...
                'EdgeColor','w','Parent',a);
            view([0,0,1]);
            c = colorbar;
            c.Color = 'w';
            c.Ticks = linspace(min(obj.A(:,1)),max(obj.A(:,1)),10);
            c.Limits = [min(obj.A(:,1)),max(obj.A(:,1))];
            c.FontSize = 20;
            colormap jet;
            title('Magnetic Vector Potential Amplitude [wb/m]',...
                'color','k','fontsize',15);
            zoom on;
        end
        function plotAmag4(obj)
            f = figure('Name','Amag','NumberTitle','off',...
                'NextPlot','new','WindowStyle','modal','Renderer','opengl',...
                'Color','k');
            a = axes(f);
            axis off equal
            patch('faces',...
                [obj.m.cl(:,[1,4,6]);obj.m.cl(:,[4,2,5]);obj.m.cl(:,[6,5,3]);obj.m.cl(:,[4,5,6])],...
                'vertices',obj.m.nodes,...
                'FaceVertexCData',obj.A,'FaceColor','interp',...
                'EdgeColor','none','Parent',a);
            view([0,0,1]);
            c = colorbar;
            c.Color = 'w';
            c.Ticks = linspace(min(obj.A(:,1)),max(obj.A(:,1)),10);
            c.Limits = [min(obj.A(:,1)),max(obj.A(:,1))];
            c.FontSize = 20;
            colormap jet;
            title('Magnetic Vector Potential Amplitude [wb/m]',...
                'color','k','fontsize',15);
            zoom on;
            hold on
            index = obj.m.edges(:,3) - obj.m.edges(:,4);
            patch('faces',obj.m.edges(logical(abs(index)),1:2),'vertices',obj.m.nodes,...
                'EdgeColor','k','Parent',a);
            set(f,'Renderer','opengl');
            zoom on;
        end
        function plotAvec(obj,scale)
            if nargin<2
                scale = 1;
            end
            close all
            quiver3(obj.m.nodes(:,1),obj.m.nodes(:,2),zeros(obj.m.Nn,1),...
                zeros(obj.m.Nn,1),zeros(obj.m.Nn,1),...
                obj.A,scale,'color','c');
            axis off equal
            view([1,1,1]);
            set(gcf,'Color','k')
            title('Magnetic Vector Potential [wb/m]',...
                'color','w','fontsize',15);
            set(gcf,'Renderer','opengl');
            shading interp
            zoom on;
        end
    end
    methods (Access = private)
        function y = getDepth(obj)
            y = obj.depth*obj.units.K_length;
        end
        function makeFalse_isElementDataAssigned(obj)
            obj.isElementDataAssigned = false;
            obj.is_Be3in_Evaluated = false;
            obj.is_Be3out_Evaluated = false;
        end
    end
end
