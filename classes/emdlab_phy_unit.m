classdef emdlab_phy_unit < handle
    
    properties
        
        unitName
        unitScaler
        isSI
        
    end
    
    methods
        
        function obj = emdlab_phy_unit(unitName, unitScaler, isSI)
            
            obj.unitName = unitName;
            obj.unitScaler = unitScaler;
            obj.isSI = isSI;
            
        end
        
    end
    
end
