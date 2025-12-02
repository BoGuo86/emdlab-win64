% EMDLAB: Electrical Machines Design Laboratorys
classdef emdlab_phy_quantity < handle
    
    properties (SetAccess = private)
        
        valueUnit
        validUnits = cell(1,[]);
        
    end
    
    methods
        
        function obj = emdlab_phy_quantity(SIunitName)
            
            obj.validUnits{end+1} = emdlab_phy_unit(SIunitName, 1, true);
            obj.valueUnit = SIunitName;
            
        end
        
        function addUnit(obj,unitName, unitScaler)
            
            if obj.isValidUnitName(unitName)
                error('Current unit name already exist.');
            end
            obj.validUnits{end+1} = emdlab_phy_unit(unitName, unitScaler, false);
            
        end
        
        function [Existence, varargout] = isValidUnitName(obj, unitName)
            
            for i = 1:numel(obj.validUnits)
                if strcmpi(obj.validUnits{i}.unitName, unitName)
                    Existence = true;
                    if nargout == 2
                        varargout{1} = i;
                    end
                    return;
                end
            end
            Existence = false;
            
        end
        
        function y = getUnitScaler(obj, unitName)
            
            [Exist, Index] = obj.isValidUnitName(unitName);
            if ~Exist
                error('Invalid unit name.');
            else
                y = obj.validUnits{Index}.unitScaler;
            end
            
        end
        
        function setValueUnit(obj, unitName)
            
            if obj.isValidUnitName(unitName)
                obj.valueUnit = unitName;
            else
                error('Invalid unit name.');
            end
            
        end
        
        function y = getValueUnit(obj)
            y = obj.valueUnit;
        end
        
        function y = getValueScaler(obj)
            y = obj.getUnitScaler(obj.valueUnit);
        end
        
        function setValue(obj, qValue)
            obj.value = qValue;
        end
        
    end
    
end
