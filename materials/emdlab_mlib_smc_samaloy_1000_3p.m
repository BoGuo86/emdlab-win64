classdef emdlab_mlib_smc_samaloy_1000_3p < handle & emdlab_mlib_electricalSteel


    methods

        function obj = emdlab_mlib_smc_samaloy_1000_3p()

            % Iron material properties (pure iron at room temperature)
            obj.ThermalConductivity.value = 30;           % W/(m·K)
            obj.HeatCapacity.value = 450;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = (1/55e-6);       % S/m

            obj.MassDensity.value = 7420;                 % kg/m³
            obj.YoungModulus.value = 170e9;              % Pa (≈ 200 GPa)
            obj.PoissonRatio.value = 0.23;                % dimensionless

            obj.gradeName = 'samaloy_1000_3p';

            obj.hb_curve = [0 	0
                46 	0.02
                75 	0.040000000000000001
                131 	0.089999999999999997
                179 	0.16
                201 	0.19
                467 	0.56000000000000005
                666 	0.79000000000000004
                1024 	1
                1801 	1.2
                4050 	1.3899999999999999
                10232 	1.5800000000000001
                24627 	1.79
                49681 	1.9399999999999999
                74681 	2.02
                99681 	2.0800000000000001
                124681 	2.1299999999999999
                149681 	2.1800000000000002
                189681 	2.25
                229681 	2.3100000000000001
                279681 	2.3900000000000001
                304681 	2.4300000000000002
                ];

            obj.evalHBCurveRelatedQuantities;

        end

    end

end