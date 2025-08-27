classdef emdlab_mlib_air < emdlab_phy_material & handle

    methods

        function obj = emdlab_mlib_air()

            % Air material properties
            obj.ThermalConductivity.value = 0.026;         % W/(m·K)
            obj.HeatCapacity.value = 1005;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;     % F/m (≈ vacuum)
            obj.ElectricConductivity.value = 0;        % S/m (essentially insulating)

            obj.MagneticPermeability.value = 1.0*4*pi*1e-7; % H/m (μ_r ≈ 1)
            obj.MagneticPermeability.isLinear = true;
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;

            obj.MassDensity.value = 1.225;                 % kg/m³
            obj.YoungModulus.value = 0;                    % Pa (air is a fluid, Young’s modulus not used)
            obj.PoissonRatio.value = 0;                    % dimensionless (not applicable for gas)

        end

    end

end