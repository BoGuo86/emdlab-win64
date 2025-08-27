classdef emdlab_mzs_common < handle
  properties
    mzs
  end
  methods
    %% Editing Functions
    function obj = emdlab_mzs_common()
    end
    function y = getDefaultMeshZoneName(obj)
      index = 0;
      mzNames = fieldnames(obj.mzs);
      while true
        index = index + 1;
        y = ['Zone',num2str(index)];
        if ~ismember(y, mzNames), break; end
      end
    end
    function mzname = checkMeshZoneExistence(obj,mzname)
      mzname = rmspaces(mzname);
      if ~isfield(obj.mzs,mzname)
        error('Specified mesh zone does not exist.');
      end
    end
    function mzname = checkMeshZoneNonExistence(obj,mzname)
      mzname = rmspaces(mzname);
      if isfield(obj.mzs,mzname)
        error('Specified mesh zone already exist.');
      end
    end
    function setMeshZoneColor(obj,mzname,color)
      mzname = obj.checkMeshZoneExistence(mzname);
      obj.mzs.(rmspaces(mzname)).color = color;
    end
    function setmzc(varargin)
      setMeshZoneColor(varargin{:});
    end
    function addmz(obj,varargin)
      if nargin<2,obj.EL.printError(0);
      elseif nargin == 2
        if ~isa(varargin{1}, 'TMZPC'),obj.EL.printError(2);end
        mzname = obj.getDefaultMeshZoneName;
        mzvalue = varargin{1};
      elseif nargin == 3
        if ~isa(varargin{1}, 'char'),obj.EL.printCharVarError('Mesh zone name');end
        mzname = obj.checkMeshZoneNonExistence(varargin{1});
        if ~isa(varargin{2}, 'TMZPC'),obj.EL.printError(2);end
        mzvalue = varargin{2};
      else,obj.EL.printError(1);
      end
      % adding new mesh zone
      obj.mzs.(mzname) = mzvalue;
      obj.mzs.(mzname).material = 'air';
      obj.mzs.(mzname).color = rand(1,3);
      % changing states
      obj.makeFalse_isGlobalMeshGenerated;
    end
    function removemz(obj,mzName)
      mzName = obj.checkMeshZoneExistence(mzName);
      delete(obj.mzs.(mzName));
      obj.mzs = rmfield(obj.mzs,mzName);
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
  end
end
