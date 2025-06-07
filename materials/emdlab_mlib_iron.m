classdef emdlab_mlib_iron < emdlab_phy_material & handle
    
    methods
        
        function obj = emdlab_mlib_iron()
            
            obj.ThermalConductivity.value = 1;
            obj.HeatCapacity.value = 1;
            obj.ElectricPermitivity.value = 0;
            obj.ElectricConductivity.value = 0;        
            obj.MagneticPermeability.value = 4000*4*pi*1e-7;
            obj.MassDensity.value = 1;
            obj.YoungModulus.value = 10;
            obj.PoissonRatio.value = 0.2;
            
        end
        
    end
    
end