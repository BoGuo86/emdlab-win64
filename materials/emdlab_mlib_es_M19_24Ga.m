classdef emdlab_mlib_es_M19_24Ga < handle & emdlab_mlib_electricalSteel


    methods

        function obj = emdlab_mlib_es_M19_24Ga()

            % Iron material properties (pure iron at room temperature)
            obj.ThermalConductivity.value = 80;           % W/(m·K)
            obj.HeatCapacity.value = 450;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = 1030000;       % S/m

            obj.MassDensity.value = 7650; % kg/m³
            obj.YoungModulus.value = 2.0e11;              % Pa (≈ 200 GPa)
            obj.PoissonRatio.value = 0.29;                % dimensionless

            obj.gradeName = 'M330-35A';

            obj.hb_curve = [0 	0
                22.280000000000001 	0.050000000000000003
                25.460000000000001 	0.10000000000000001
                31.829999999999998 	0.14999999999999999
                47.740000000000002 	0.35999999999999999
                63.659999999999997 	0.54000000000000004
                79.569999999999993 	0.65000000000000002
                159.15000000000001 	0.98999999999999999
                318.30000000000001 	1.2
                477.45999999999998 	1.28
                636.61000000000001 	1.3300000000000001
                795.76999999999998 	1.3600000000000001
                1591.5 	1.4399999999999999
                3183 	1.52
                4774.6000000000004 	1.5800000000000001
                6366.1000000000004 	1.6299999999999999
                7957.6999999999998 	1.6699999999999999
                15915 	1.8
                31830 	1.8999999999999999
                111407 	2
                190984 	2.1000000000000001
                350138 	2.2999999999999998
                509292 	2.5];

            obj.evalHBCurveRelatedQuantities;

        end

    end

end