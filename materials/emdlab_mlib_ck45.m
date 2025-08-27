classdef emdlab_mlib_ck45 < emdlab_phy_material

    methods

        function obj = emdlab_mlib_ck45()

            % CK45 steel material properties (medium-carbon steel, room temperature)
            obj.ThermalConductivity.value = 50;           % W/(m·K)
            obj.HeatCapacity.value = 470;                 % J/(kg·K)
            obj.ElectricPermitivity.value = 8.854e-12;    % F/m (≈ vacuum, usually ignored for conductors)
            obj.ElectricConductivity.value = 6.99e6;      % S/m

            obj.MagneticPermeability.value = 1000*4*pi*1e-7; % H/m (μ_r ≈ 1000 for low-carbon steel, nonlinear)
            obj.MagneticPermeability.isLinear = false;    % CK45 is ferromagnetic
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isScalar = true;

            obj.MassDensity.value = 7800;                 % kg/m³
            obj.YoungModulus.value = 2.1e11;              % Pa (≈ 210 GPa)
            obj.PoissonRatio.value = 0.29;                % dimensionless

        end

    end

end
