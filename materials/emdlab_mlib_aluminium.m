classdef emdlab_mlib_aluminium < emdlab_phy_material
    
    methods
        
        function obj = emdlab_mlib_aluminium()
            
            obj.ThermalConductivity.value = 380;
            obj.HeatCapacity.value = 400;
            obj.ElectricPermitivity.value = 0;
            obj.ElectricConductivity.value = 0;
            obj.MagneticPermeability.value = 0.999*4*pi*1e-7;
            obj.MassDensity.value = 8500;
            obj.YoungModulus.value = 10;
            obj.PoissonRatio.value = 0.2;
            
        end
        
    end
    
end