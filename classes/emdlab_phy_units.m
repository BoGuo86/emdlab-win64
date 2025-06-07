classdef emdlab_phy_units < handle
    
    properties (SetAccess = private)
        
        length
        time
        current
        voltage
        currentDensity
        magneticVectorPotential
        magnetisation
        frequency
        surfaceLossDensity
        volumeLossDensity
        temperature
        convectionHeatTransferCoefficient
        power
        surfaceForceDensity
        
    end
    
    properties (Dependent = true)
        
        k_length
        k_time
        k_current
        k_voltage
        k_currentDensity
        k_magneticVectorPotential
        k_frequency
        k_surfaceLossDensity
        k_volumeLossDensity
        k_temperature
        k_convectionHeatTransfer
        k_magnetisation
        k_power
        k_convectionHeatTransferCoefficient
        k_surfaceForceDensity
        
    end
    
    methods
        
        % constructor
        function obj = emdlab_phy_units()
            
            % length
            obj.length = emdlab_phy_quantity('m');
            obj.length.addUnit('um', 1e-6);
            obj.length.addUnit('mm', 1e-3);
            obj.length.addUnit('km', 1e3);
            obj.length.addUnit('cm', 1e-2);
            
            % time
            obj.time = emdlab_phy_quantity('s');
            obj.time.addUnit('ns', 1e-9);
            obj.time.addUnit('us', 1e-6);
            obj.time.addUnit('ms', 1e-3);
            
            % current
            obj.current = emdlab_phy_quantity('a');
            obj.current.addUnit('ua', 1e-6);
            obj.current.addUnit('ma', 1e-3);
            obj.current.addUnit('ka', 1e3);
            
            % voltage
            obj.voltage = emdlab_phy_quantity('v');
            obj.voltage.addUnit('uv', 1e-6);
            obj.voltage.addUnit('mv', 1e-3);
            obj.voltage.addUnit('kv', 1e3);
            
            % currentDensity
            obj.currentDensity = emdlab_phy_quantity('a/m^2');
            obj.currentDensity.addUnit('a/mm^2', 1e6);
            
            % magneticVectorPotential
            obj.magneticVectorPotential = emdlab_phy_quantity('a/m');
            obj.magneticVectorPotential.addUnit('a/mm', 1e3);
            
            % magneticVectorPotential
            obj.frequency = emdlab_phy_quantity('hz');
            obj.frequency.addUnit('mhz', 1e-3);
            obj.frequency.addUnit('khz', 1e3);
            
            % surfaceLossDensity
            obj.surfaceLossDensity = emdlab_phy_quantity('w/m^2');
            obj.surfaceLossDensity.addUnit('w/mm^2', 1e6);
            
            % volumeLossDensity
            obj.volumeLossDensity = emdlab_phy_quantity('w/m^3');
            obj.volumeLossDensity.addUnit('w/mm^3', 1e9);
            
            % temperature
            obj.temperature = emdlab_phy_quantity('c');
            
            % convection heat transfer
            obj.convectionHeatTransferCoefficient = emdlab_phy_quantity('w/(m^2*c)');
            obj.convectionHeatTransferCoefficient.addUnit('w/(mm^2*c)', 1e6);
            
            % magnetisation
            obj.magnetisation = emdlab_phy_quantity('a/m');
            obj.magnetisation.addUnit('a/mm', 1e3);
            
            % power
            obj.power = emdlab_phy_quantity('w');
            obj.power.addUnit('mw', 1e-3);
            obj.power.addUnit('kw', 1e3);
            
            % surfaceForceDensity
            obj.surfaceForceDensity = emdlab_phy_quantity('n/m^2');
            obj.surfaceForceDensity.addUnit('n/mm^2', 1e6);
            
        end
        
        % operations
        function y = isQuantity(obj, quantityName)
            quantityName = strrep(quantityName, ' ', '');
            qNames = fieldnames(obj);
            for i = 1:numel(qNames)
                if strcmpi(quantityName, qNames{i})
                    y = true;
                    return;
                end
            end
            y = false;
        end
        
        function y = getUnitScaler(obj, quantityName, quantityUnit)
            if ~obj.isQuantity(quantityName)
                error('Invalid Quantity name.');
            end
            y = obj.(quantityName).getUnitScaler(quantityUnit);
        end
        
        function setQuantityUnit(obj, quantityName, quantityUnit)
            if ~obj.isQuantity(quantityName)
                error('Invalid Quantity name.');
            end
            quantityUnit = strrep(quantityUnit, ' ', '');
            quantityUnit = lower(quantityUnit);
            obj.(quantityName).setValueUnit(quantityUnit);
        end
        
        function y = getQuantityScaler(obj, quantityName)
            if ~obj.isQuantity(quantityName)
                error('Invalid Quantity name.');
            end
            y = obj.(quantityName).getValueScaler();
        end
        
        % getters
        function y = get.k_length(obj)
            y = obj.getQuantityScaler('length');
        end
        function y = get.k_time(obj)
            y = obj.getQuantityScaler('time');
        end
        function y = get.k_current(obj)
            y = obj.getQuantityScaler('current');
        end
        function y = get.k_currentDensity(obj)
            y = obj.getQuantityScaler('currentDensity');
        end
        function y = get.k_magneticVectorPotential(obj)
            y = obj.getQuantityScaler('magneticVectorPotential');
        end
        function y = get.k_frequency(obj)
            y = obj.getQuantityScaler('frequency');
        end
        function y = get.k_surfaceLossDensity(obj)
            y = obj.getQuantityScaler('surfaceLossDensity');
        end
        function y = get.k_volumeLossDensity(obj)
            y = obj.getQuantityScaler('volumeLossDensity');
        end
        function y = get.k_temperature(obj)
            y = obj.getQuantityScaler('temperature');
        end
        function y = get.k_voltage(obj)
            y = obj.getQuantityScaler('voltage');
        end
        function y = get.k_magnetisation(obj)
            y = obj.getQuantityScaler('magnetisation');
        end
        function y = get.k_power(obj)
            y = obj.getQuantityScaler('power');
        end
        function y = get.k_convectionHeatTransferCoefficient(obj)
            y = obj.getQuantityScaler('convectionHeatTransferCoefficient');
        end
        function y = get.k_surfaceForceDensity(obj)
            y = obj.getQuantityScaler('surfaceForceDensity');
        end
    end
end
