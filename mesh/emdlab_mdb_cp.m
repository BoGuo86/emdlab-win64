classdef emdlab_mdb_cp < handle
    % common properties of mesh databases
    
    properties
        
        % mesh zones
        mzs (1,1) struct;
        
        % materials
        mts (1,1) struct;
        
        % states
        isGlobalMeshGenerated (1,1) logical = false;
        
    end
    
    properties (Dependent = true)
        
        % number of mesh zones
        Nmzs (1,1) double {mustBeNonnegative, mustBeInteger} = 0;
        
    end
    
    methods
        
        function obj = emdlab_mdb_cp()
        end
        
        function y = get.Nmzs(obj)
            
            y = numel(fieldnames(obj.mzs));
            
        end
        
        function y = getDefaultMeshZoneName(obj)
            
            index = 0;
            mzNames = fieldnames(obj.mzs);
            
            while true
                index = index + 1;
                y = ['Zone', num2str(index)];
                if ~ismember(y, mzNames), break; end
            end
            
        end
        
        function mzName = checkMeshZoneExistence(obj, mzName)
            
            mzName = rmspaces(mzName);
            
            if ~isfield(obj.mzs, mzName)
                error('Specified mesh zone does not exist.');
            end
            
        end
        
        function mzName = checkMeshZoneNonExistence(obj, mzName)
            
            mzName = rmspaces(mzName);
            
            if isfield(obj.mzs, mzName)
                error('Specified mesh zone already exist.');
            end
            
        end
        
        function setMeshZoneColor(obj, mzName, R, G, B)
            
            mzName = obj.checkMeshZoneExistence(mzName);
            obj.mzs.(rmspaces(mzName)).color = [R,G,B]/255;
            
        end
        
        function setmzc(varargin)
            
            setMeshZoneColor(varargin{:});
            
        end
        
        function addMeshZone(obj, varargin)
            
            if nargin == 2
                
                if ~isa(varargin{1}, 'emdlab_m2d_tmz')
                    error('Mesh zone class must be <emdlab_m2d_tmz>.');
                end
                
                mzName = obj.getDefaultMeshZoneName;
                mzptr = varargin{1};
                
            elseif nargin == 3
                mzName = obj.checkMeshZoneNonExistence(varargin{1});
                
                if ~isa(varargin{2}, 'emdlab_m2d_tmz')
                    error('Mesh zone class must be <emdlab_m2d_tmz>.');
                end
                
                mzptr = varargin{2};
            else
                error('Wrong number of arguments.');
            end
            
            % adding new mesh zone
            obj.mzs.(mzName) = mzptr;
            obj.mzs.(mzName).material = 'air';
            obj.mzs.(mzName).color = rand(1, 3);
            
            % changing states
            obj.makeFalse_isGlobalMeshGenerated;
            
        end
        
        function addmz(obj, varargin)
            
            obj.addMeshZone(varargin{:});
            
        end
        
        function removeMeshZone(obj, mzName)
            
            mzName = obj.checkMeshZoneExistence(mzName);
            % delete(obj.mzs.(mzName));
            obj.mzs = rmfield(obj.mzs, mzName);
            % changing states
            obj.makeFalse_isGlobalMeshGenerated;
            
        end
        
        function removemz(obj, varargin)
            
            obj.removeMeshZone(varargin{:});
            
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
        
        %% Material Library
        
        function setMaterial(obj, mzName, mName)
            
            mzName = obj.checkMeshZoneExistence(mzName);
            mName = strrep(mName, ' ', '');
            
            if ~ismember(mName, fieldnames(obj.mts))
                error(['Material <<', mName, '>> does not found.']);
            end
            
            obj.mzs.(mzName).material = mName;
            
        end
        
        function addMaterial(obj, mName, mObject)
            
            if nargin == 2
                obj.mts.(mName) = eval(['emdlab_mlib_', mName]);
            else
                obj.mts.(mName) = mObject;
            end
            
        end
        
    end
    
end