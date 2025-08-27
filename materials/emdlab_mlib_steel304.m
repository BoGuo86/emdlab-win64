classdef emdlab_mlib_steel304 < emdlab_phy_material

    methods

        function obj = emdlab_mlib_steel304()

            % Stainless Steel 304 material properties (at room temperature)
            obj.ThermalConductivity.value = 16;           % W/(m·K)
            obj.HeatCapacity.value = 500;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = 1.45e6;      % S/m (depends on alloy and treatment)

            obj.MagneticPermeability.value = 1.0*4*pi*1e-7; % H/m (μ_r ≈ 1, 304 is generally non-magnetic)
            obj.MagneticPermeability.isLinear = true;
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;

            obj.MassDensity.value = 8000;                 % kg/m³
            obj.YoungModulus.value = 1.9e11;              % Pa (≈ 190 GPa)
            obj.PoissonRatio.value = 0.30;                % dimensionless

        end

    end

end
