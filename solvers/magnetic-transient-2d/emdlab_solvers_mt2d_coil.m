% two-dimensional magnetic-transient coil
% a coil is made of a series of coil arms: there is no parallel branches

classdef emdlab_solvers_mt2d_coil < handle & matlab.mixin.SetGet

    properties

        % the number of coil turns
        coilArms (1,:) string;

        % coil index
        ci (1,1) double;

        % coil arm directions: 'positive' -> z+, or 'negative' -> z-
        directions (1,:) double;

        % voltage
        voltage (1,:) double = 0;

        % current
        current (1,:) double = 0;

        % induced voltage
        inducedVoltage (1,:) double = 0;

        % flux linkage
        fluxLinkage (1,:) double = 0;

        % dc resistance of the coil
        Rdc (1,1) double;      

        % initial current of coil, is ued in voltage fed case (or external)
        initialCurrent (1,1) double; 

        % fed type is 'voltage' or 'current'
        fedType (1,:) char = 'voltage';

        % eddy type is 'stranded' or 'solid'
        eddyType (1,:) char = 'stranded';

        % external resistance
        Rex (1,1) double {mustBeNonnegative} = 0;

        % external inductance
        Lex (1,1) double {mustBeNonnegative} = 0;

        % Qvec to calculate coil flux linkage
        Qvec (1,:) double;

        % voltage or current value
        vcValue (1,1) = 0;

        % flags
        isCageMember = false;        
        isStarConnectionMember = false;

    end

    properties (Dependent = true)

        % number of coil arms
        NcoilArms (1,1) double;

    end

    methods

        function obj = emdlab_solvers_mt2d_coil()
        end   

        function set.eddyType(obj, newValue)
            if ~ismember(newValue, ["stranded", "solid"])
                error('The coil eddy type must be <solid> or <stranded>.');
            end
            obj.eddyType = newValue;
        end

        function set.fedType(obj, newValue)
            if ~ismember(newValue, ["voltage", "current"])
                error('The coil source type must be <voltage> or <current>.');
            end
            obj.fedType = newValue;
        end

        function addCoilArm(obj, newCoilArmName, direction)

            obj.coilArms(end+1) = newCoilArmName;
            obj.directions(end+1) = direction;

        end

        function y = get.NcoilArms(obj)
            y = numel(obj.coilArms);
        end

        function y = isCurrentFed(obj)
            if strcmpi(obj.fedType, 'current')
                y = true;
            else
                y = false;
            end
        end

        function y = isVoltageFed(obj)
            if strcmpi(obj.fedType, 'voltage')
                y = true;
            else
                y = false;
            end
        end

        function y = isStranded(obj)
            if strcmpi(obj.eddyType, 'stranded')
                y = true;
            else
                y = false;
            end
        end

        function y = isSolid(obj)
            if strcmpi(obj.eddyType, 'solid')
                y = true;
            else
                y = false;
            end
        end

        function setVoltage(obj, value)
            if obj.isCurrentFed
                error('This coil is current fed. You must set the coil current.');
            end
            if isa(value, 'function_handle') || (isscalar(value) && isnumeric(value))
                obj.vcValue = value;
            else
                error('Coil voltage must be a function handle or a scalar numeric value.');
            end
        end

        function setCurrent(obj, value)
            if obj.isVoltageFed
                error('This coil is voltage fed. You must set the coil voltage.');
            end
            if isa(value, 'function_handle') || (isscalar(value) && isnumeric(value))
                obj.vcValue = value;
            else
                error('Coil current must be a function handle or a scalar numeric value.');
            end
        end

        function setInitialCurrent(obj, value)
            if obj.isCurrentFed
                error('This coil is current fed. Initial current can be set for voltage coil.');
            end
            if isscalar(value) && isnumeric(value)
                obj.voltageValue = value;
            else
                error('Initial current of a voltage coil must be a scalar numeric value.');
            end
        end

        function y = getCurrent(obj, timePoint)
            if obj.isVoltageFed
                error('You cannot get the current of a voltage fed coil.');
            else
                if isa(obj.vcValue, 'function_handle')
                    y = obj.vcValue(timePoint);
                else
                    y = obj.vcValue;
                end
            end
            if isnan(y)
                error('The coil current is non. Check the value of coil current.');
            end
        end

        function y = getVoltage(obj, timePoint)
            if obj.isCurrentFed
                error('You cannot get the voltage of a current fed coil.');
            else
                if isa(obj.vcValue, 'function_handle')
                    y = obj.vcValue(timePoint);
                else
                    y = obj.vcValue;
                end
            end
            if isnan(y)
                error('The coil voltage is non. Check the value of coil voltage.');
            end
        end

    end
end
