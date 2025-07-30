classdef emdlab_mlib_pm_Y22H < handle & emdlab_mlib_permanentMagnet
    
    methods
        
        function obj = emdlab_mlib_pm_Y22H()
            obj.ThermalConductivity.value = 1;
            obj.HeatCapacity.value = 1;
            obj.ElectricPermitivity.value = 0;
            obj.ElectricConductivity.value = 0;        
            obj.MagneticPermeability.value = 4*pi*1e-7;
            obj.MassDensity.value = 1;
            obj.YoungModulus.value = 10;
            obj.PoissonRatio.value = 0.2;
        end
        
    end
    
end