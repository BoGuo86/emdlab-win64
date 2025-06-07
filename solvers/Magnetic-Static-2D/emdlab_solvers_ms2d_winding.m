% developer: https://ComProgExpert.com
% winding object in magnetic-static analysis

classdef emdlab_solvers_ms2d_winding < handle & matlab.mixin.SetGet

    properties

        % current amplitude
        current(1, 1) double = 0;
        
        % number of parallel path
        np(1, 1) double {mustBePositive, mustBeInteger} = 1;

    end

    properties (SetAccess = private)

        % mesh zones assigned to the winding
        mzsName cell = {};

    end

    properties (Dependent = true)

        % number of matrix mesh zones
        Nmzs(1, 1) double;

    end

    methods

        function obj = emdlab_solvers_ms2d_winding(varargin)

            if numel(varargin)
                set(obj, varargin{:});
            end
            
        end

        function addMeshZone(obj, mzName)

            % chekers
            mzName = obj.checkMeshZoneNonExistence(mzName);
            obj.mzsName{end + 1} = mzName;

        end

        function mzName = checkMeshZoneExistence(obj, mzName)

            if ~ismember(mzName, obj.mzsName)
                error('Specified mesh zone does not exist.');
            end

        end

        function mzName = checkMeshZoneNonExistence(obj, mzName)

            if ismember(mzName, obj.mzsName)
                error('Specified mesh zone already exist.');
            end

        end

        function y = get.Nmzs(obj)
            y = numel(obj.mzsName);
        end

    end

end
