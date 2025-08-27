classdef emdlab_mlib_zinc < emdlab_phy_material

    methods

        function obj = emdlab_mlib_zinc()

            % Zinc material properties (at room temperature)
            obj.ThermalConductivity.value = 116;          % W/(m·K)
            obj.HeatCapacity.value = 390;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = 1.69e7;      % S/m

            obj.MagneticPermeability.value = 1.0*4*pi*1e-7; % H/m (μ_r ≈ 1)
            obj.MagneticPermeability.isLinear = true;
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;

            obj.MassDensity.value = 7135;                 % kg/m³
            obj.YoungModulus.value = 1.0e11;              % Pa (≈ 100 GPa)
            obj.PoissonRatio.value = 0.25;                % dimensionless

        end

    end

end
