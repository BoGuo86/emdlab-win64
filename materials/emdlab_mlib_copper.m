classdef emdlab_mlib_copper < emdlab_phy_material
    
    methods
        
        function obj = emdlab_mlib_copper()
            
            obj.ThermalConductivity.value = 380;
            obj.HeatCapacity.value = 400;
            obj.ElectricPermitivity.value = 0;
            obj.ElectricConductivity.value = 0;
            obj.MagneticPermeability.value = 0.999*4*pi*1e-7;
            obj.MagneticPermeability.isLinear = true;
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;
            
            obj.MassDensity.value = 8500;
            obj.YoungModulus.value = 10;
            obj.PoissonRatio.value = 0.2;
            
        end
        
    end
    
end