classdef emdlab_mlib_es_steel1008 < handle & emdlab_mlib_electricalSteel
    
    
    methods
        
        function obj = emdlab_mlib_es_steel1008()
            
            % Iron material properties (pure iron at room temperature)
            obj.ThermalConductivity.value = 80;           % W/(m·K)
            obj.HeatCapacity.value = 450;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = 1030000;       % S/m

            obj.MassDensity.value = 7870;                 % kg/m³
            obj.YoungModulus.value = 2.0e11;              % Pa (≈ 200 GPa)
            obj.PoissonRatio.value = 0.29;                % dimensionless
            
            obj.gradeName = 'Steel-1008';
            
            obj.hb_curve = [0 0
                159.19999999999999 	0.2402
                318.30000000000001 	0.86539999999999995
                477.5 	1.1106
                636.60000000000002 	1.2458
                795.79999999999995 	1.331
                1591.5 	1.5
                3183.0999999999999 	1.6000000000000001
                4774.6000000000004 	1.6830000000000001
                6366.1999999999998 	1.7410000000000001
                7957.6999999999998 	1.78
                15915.5 	1.905
                31831 	2.0249999999999999
                47746.5 	2.085
                63662 	2.1299999999999999
                79577.5 	2.165
                159155 	2.2799999999999998
                318310 	2.4849999999999999
                397887 	2.5851000000000002];
            
        end
        
    end
    
end