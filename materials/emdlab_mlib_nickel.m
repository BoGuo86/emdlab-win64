classdef emdlab_mlib_nickel < emdlab_phy_material

    methods

        function obj = emdlab_mlib_nickel()

            % Nickel material properties (pure nickel at room temperature)
            obj.ThermalConductivity.value = 90;           % W/(m·K)
            obj.HeatCapacity.value = 440;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = 1.43e7;      % S/m

            obj.MagneticPermeability.value = 600*4*pi*1e-7; % H/m (relative permeability μ_r ≈ 600 for nickel)
            obj.MagneticPermeability.isLinear = false;    % Nickel is nonlinear in magnetic response
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;

            obj.MassDensity.value = 8900;                 % kg/m³
            obj.YoungModulus.value = 2.0e11;              % Pa (≈ 200 GPa)
            obj.PoissonRatio.value = 0.31;                % dimensionless

        end

    end

end
