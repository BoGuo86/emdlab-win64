classdef emdlab_mlib_steel1010 < emdlab_phy_material

    methods

        function obj = emdlab_mlib_steel1010()

            % Steel 1010 material properties (low-carbon steel, room temperature)
            obj.ThermalConductivity.value = 50;           % W/(m·K)
            obj.HeatCapacity.value = 460;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = 6.99e6;      % S/m

            obj.MagneticPermeability.value = 1000*4*pi*1e-7; % H/m (μ_r ≈ 1000, ferromagnetic)
            obj.MagneticPermeability.isLinear = false;    % Nonlinear magnetic response
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;

            obj.MassDensity.value = 7850;                 % kg/m³
            obj.YoungModulus.value = 2.0e11;              % Pa (≈ 200 GPa)
            obj.PoissonRatio.value = 0.29;                % dimensionless

        end

    end

end
