classdef emdlab_solvers_IHLMSTL3 < handle & emdlab_MSTL3
  methods
    %% Initialization
    function obj = emdlab_solvers_IHLMSTL3(m)
      obj.m = m;
      % default properties of mzs
      mzNames = fieldnames(obj.m.mzs);
      for i = 1:numel(mzNames)
        obj.setdp(mzNames{i});
      end
    end
    %% Solver
    function assignEdata(obj)
      if obj.isElementDataAssigned, return; end
      % preparing mesh data
      obj.m.evalKexy1_TL3;
      tic, disp('-------------------------------------------------------');
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
                obj.m.mzs.(mzname).props.excitation.value/obj.m.mzs.(mzname).getArea;
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
    function obj = solve(obj)
      obj.assignEdata;
      % Construction of [K] and [F]
      tic
      disp('*******************************************************')
      % Assembeling [F]
      fi = repmat(obj.edata.InternalCurrentDensity,3,1);
      Mx = repmat(obj.edata.Magnetization(1,:),3,1)*sparse(1:obj.m.Ne,1:obj.m.Ne,obj.m.gea);
      My = repmat(obj.edata.Magnetization(2,:),3,1)*sparse(1:obj.m.Ne,1:obj.m.Ne,obj.m.gea);
      F = (fi.*obj.m.Fe*obj.units.K_currentDensity+...
        (obj.m.gphiy.*Mx-obj.m.gphix.*My)/obj.units.K_length);
      % applying scales on load vector
      F = F*obj.units.K_length^2/obj.units.K_magneticVectorPotential;
      F = sparse(obj.m.cl',ones(3*obj.m.Ne,1),F);
      % Assembeling [K]
      [Iindex,Jindex] = getij(3,1);
      K = sparse(obj.m.cl(:,Iindex)',...
        obj.m.cl(:,Jindex)',...
        repmat(obj.edata.MagneticReluctivity,9,1).* ...
        obj.m.Ke(getkindex(3),:));
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
      obj.evalB;
      disp('initial geuss evaluated.')
      toc
    end
  end
end
