classdef emdlab_mlib_brass < emdlab_phy_material

    methods

        function obj = emdlab_mlib_brass()

            % Brass material properties (typical values at room temperature)
            obj.ThermalConductivity.value = 109;          % W/(m·K)
            obj.HeatCapacity.value = 380;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = 1.5e7;       % S/m (depends on alloy composition)

            obj.MagneticPermeability.value = 1.0*4*pi*1e-7; % H/m (μ_r ≈ 1)
            obj.MagneticPermeability.isLinear = true;
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;

            obj.MassDensity.value = 8500;                 % kg/m³ (typical for common brass alloys)
            obj.YoungModulus.value = 1.0e11;              % Pa (≈ 100 GPa)
            obj.PoissonRatio.value = 0.34;                % dimensionless

        end

    end

end
