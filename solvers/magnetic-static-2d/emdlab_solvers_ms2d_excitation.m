% EMDLAB: Electrical Machines Design Laboratory

classdef emdlab_solvers_ms2d_excitation < handle & matlab.mixin.SetGet

    properties

        % type of magnetic-static excitation, (1) current, (2) currentDensity
        type (1,:) char;

        % the value of excitation
        value (1,1) double = 0;

    end

    methods

        function obj = emdlab_solvers_ms2d_excitation(varargin)

            % default types
            obj.value = 1;
            obj.type = 'current';
            set(obj, varargin{:});

        end

        function set.type(obj, xtype)
            if ischar(xtype)
                xtype = rmspaces(xtype);
                if strcmpi(xtype,'c') || strcmpi(xtype,'current')
                    obj.type = 'current';
                elseif strcmpi(xtype,'cd') || strcmpi(xtype,'currentDensity')
                    obj.type = 'currentDensity';
                else
                    error('Excitation type must be <<current[or c]>> or <<currentDensity [or cd]>>.');
                end
            else
                error('Excitation type must be a string.');
            end
        end

    end

end
