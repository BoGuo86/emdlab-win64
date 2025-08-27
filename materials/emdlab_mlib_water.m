classdef emdlab_mlib_water < emdlab_phy_material

    methods

        function obj = emdlab_mlib_water()

            % Water material properties (at room temperature)
            obj.ThermalConductivity.value = 0.6;           % W/(m·K)
            obj.HeatCapacity.value = 4186;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 7.08e-10;      % F/m (≈ ε0 * 80)
            obj.ElectricConductivity.value = 5.5e-6;       % S/m (pure water very low conductivity)

            obj.MagneticPermeability.value = 1.0*4*pi*1e-7; % H/m (μ_r ≈ 1)
            obj.MagneticPermeability.isLinear = true;
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;

            obj.MassDensity.value = 1000;                  % kg/m³
            obj.YoungModulus.value = 2.2e9;                % Pa (≈ 2.2 GPa, bulk modulus equivalent for FEM)
            obj.PoissonRatio.value = 0.499;                % dimensionless (nearly incompressible)

        end

    end

end
