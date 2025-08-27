classdef emdlab_phy_materialProperty
    
    properties
        
        % value of material property: a two dimensional array
        value (1,1) double;
        
        % scalar vs vector
        isScalar (1,1) logical;
        
        % linear vs non-linear
        isLinear (1,1) logical;
        
        % isotropic vs non-isotropic
        isIsotropic (1,1) logical;
        
    end
    
    methods
        
        function obj = emdlab_phy_materialProperty()
            
            obj.value = 0;
            obj.isLinear = true;
            obj.isScalar = true;
            obj.isIsotropic = true;
            
        end
        
    end
end
