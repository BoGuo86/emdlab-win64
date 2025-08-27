classdef emdlab_mlib_steel316 < emdlab_phy_material

    methods

        function obj = emdlab_mlib_steel316()

            % Stainless Steel 316 material properties (at room temperature)
            obj.ThermalConductivity.value = 16;           % W/(m·K)
            obj.HeatCapacity.value = 500;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = 1.0e6;       % S/m (depends on alloy composition)

            obj.MagneticPermeability.value = 1.0*4*pi*1e-7; % H/m (μ_r ≈ 1, 316 is non-magnetic)
            obj.MagneticPermeability.isLinear = true;
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;

            obj.MassDensity.value = 8000;                 % kg/m³
            obj.YoungModulus.value = 2.0e11;              % Pa (≈ 200 GPa)
            obj.PoissonRatio.value = 0.30;                % dimensionless

        end

    end

end
