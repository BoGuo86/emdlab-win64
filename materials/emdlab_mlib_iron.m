classdef emdlab_mlib_iron < emdlab_phy_material

    methods

        function obj = emdlab_mlib_iron()

            % Iron material properties (pure iron at room temperature)
            obj.ThermalConductivity.value = 80;           % W/(m·K)
            obj.HeatCapacity.value = 450;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = 1030000;       % S/m

            obj.MagneticPermeability.value = 4000*4*pi*1e-7; % H/m (relative permeability μ_r ≈ 5000 for low-carbon iron)
            obj.MagneticPermeability.isLinear = true;    % Iron is nonlinear in magnetic response
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;

            obj.MassDensity.value = 7870;                 % kg/m³
            obj.YoungModulus.value = 2.0e11;              % Pa (≈ 200 GPa)
            obj.PoissonRatio.value = 0.29;                % dimensionless

        end

    end

end
