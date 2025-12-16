classdef emdlab_mlib_es_M330_35A_AM < handle & emdlab_mlib_electricalSteel
    
    
    methods
        
        function obj = emdlab_mlib_es_M330_35A_AM()
            
            % Iron material properties (pure iron at room temperature)
            obj.ThermalConductivity.value = 80;           % W/(m·K)
            obj.HeatCapacity.value = 450;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = 1030000;       % S/m

            obj.MassDensity.value = 7870;                 % kg/m³
            obj.YoungModulus.value = 2.0e11;              % Pa (≈ 200 GPa)
            obj.PoissonRatio.value = 0.29;                % dimensionless
            
            obj.gradeName = 'M330-35A';
            
            obj.hb_curve = [0 	0
5 	0.0081679000000000005
10 	0.017769199999999999
15 	0.029215999999999999
25 	0.058117799999999997
33.399999999999999 	0.10000000000000001
43.600000000000001 	0.20000000000000001
50.799999999999997 	0.29999999999999999
57.200000000000003 	0.40000000000000002
63.600000000000001 	0.5
70.400000000000006 	0.59999999999999998
78.099999999999994 	0.69999999999999996
87.200000000000003 	0.80000000000000004
98.700000000000003 	0.90000000000000002
114 	1
136 	1.1000000000000001
172 	1.2
242 	1.3
428 	1.3999999999999999
1027 	1.5
2576 	1.6000000000000001
5409 	1.7
9677 	1.8
18943 	1.9166000000000001
37081 	1.9930000000000001
72586 	2.0649999999999999
142090 	2.1663999999999999
278140 	2.3445
544460 	2.6827999999999999
1065800 	3.3397999999999999];
            
            obj.evalHBCurveRelatedQuantities;
            
        end
        
    end
    
end