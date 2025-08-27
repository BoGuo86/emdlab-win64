classdef emdlab_mlib_ck60 < emdlab_phy_material

    methods

        function obj = emdlab_mlib_ck60()

            % CK60 steel material properties (high-carbon steel, room temperature)
            obj.ThermalConductivity.value = 50;           % W/(m·K)
            obj.HeatCapacity.value = 460;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = 6.99e6;      % S/m

            obj.MagneticPermeability.value = 1300*4*pi*1e-7; % H/m (μ_r ≈ 1300, ferromagnetic)
            obj.MagneticPermeability.isLinear = false;    % Nonlinear magnetic response
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;

            obj.MassDensity.value = 7850;                 % kg/m³
            obj.YoungModulus.value = 2.1e11;              % Pa (≈ 210 GPa)
            obj.PoissonRatio.value = 0.29;                % dimensionless

        end

    end

end
