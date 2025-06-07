%% EMDLAB
% mfile description:
% this mfile write materials in .dat binary files
% for each material at first physical properties must be set and then the
% function WriteMDataBinary() must be called.
% physical properties:
% ===> ThermalConductivity
% ===> HeatCapacity
% ===> ElectricPermitivity
% ===> ElectricConductivity
% ===> MagneticPermeability
% ===> MassDensity
% ===> YoungModulus
% ===> PoissonRatio

function WriteMaterials()

% directory that materials data will be saved
fdir = [cd,'\materials\'];

% physical constant class
pcts = emdlab_phy_constants;

%% air
m.ThermalConductivity = 1;
m.HeatCapacity = 1;
m.ElectricPermitivity = pcts.e0;
m.ElectricConductivity = 0;
m.MagneticPermeability = pcts.mu0;
m.MassDensity = 1;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'air',m)

%% copper
m.ThermalConductivity = 380;
m.HeatCapacity = 400;
m.ElectricPermitivity = 8.85*1e-12;
m.ElectricConductivity = 58e6;
m.MagneticPermeability = 0.999*pcts.mu0;
m.MassDensity = 450;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'copper',m)

%% iron
m.ThermalConductivity = 400;
m.HeatCapacity = 1;
m.ElectricPermitivity = pcts.e0;
m.ElectricConductivity = 50e5;
m.MagneticPermeability = 4000*pcts.mu0;
m.MassDensity = 7800;
m.YoungModulus = 10e9;
m.PoissonRatio = 0.23;
WriteMDataBinary(fdir,'iron',m)

%% laminatedIron
m.ThermalConductivity = 1;
m.HeatCapacity = 1;
m.ElectricPermitivity = pcts.e0;
m.ElectricConductivity = 0;
m.MagneticPermeability = 4000*pcts.mu0;
m.MassDensity = 1;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'laminatedIron',m)

%% lamination
m.ThermalConductivity = [30 30 0.6];
m.HeatCapacity = 2;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 0;
m.MagneticPermeability = 4000*pcts.mu0;
m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'lamination',m)

%% nilamination
m.ThermalConductivity = [30,30,0.6];
m.HeatCapacity = 2;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 0;
m.MagneticPermeability = 4000*pcts.mu0;
m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'nilamination',m)

%% wedge
m.ThermalConductivity = 0.3;
m.HeatCapacity = 2;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 0;
m.MagneticPermeability = 4000*pcts.mu0;
m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'wedge',m)

%% eqwin
m.ThermalConductivity = [0.83,0.83,50];
m.HeatCapacity = 2;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 0;
m.MagneticPermeability = 4000*pcts.mu0;
m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'eqwin',m)

%% eqwin
m.ThermalConductivity = 0.83;
m.HeatCapacity = 2;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 0;
m.MagneticPermeability = 4000*pcts.mu0;
m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'eqwin2d',m)

%% ew
m.ThermalConductivity = 2;
m.HeatCapacity = 2;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 0;
m.MagneticPermeability = 4000*pcts.mu0;
m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'ew',m)

%% material name: m19_24g
m.ThermalConductivity = 30;
m.HeatCapacity = 1;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 0;

f = fopen([cd,'\materials\hb-curves\','m19_24g.tab'],'r');
HB = fscanf(f,'%f');
fclose(f);
m.MagneticPermeability = [HB(1:2:end),HB(2:2:end)];

m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'m19_24g',m);

%% material name: m350_50a
m.ThermalConductivity = 30;
m.HeatCapacity = 1;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 0;

f = fopen([cd,'\materials\hb-curves\','m350_50a.tab'],'r');
HB = fscanf(f,'%f');
fclose(f);
m.MagneticPermeability = [HB(1:2:end),HB(2:2:end)];

m.MassDensity = 7850;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'m350_50a',m);


%% material name: M530-50A
m.ThermalConductivity = 30;
m.HeatCapacity = 1;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 0;

f = fopen([cd,'\materials\hb-curves\','m530_50a.tab'],'r');
HB = fscanf(f,'%f');
fclose(f);
m.MagneticPermeability = [HB(1:2:end),HB(2:2:end)];

m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'m530_50a',m)

%% ===>  steel_1008
m.ThermalConductivity = 30;
m.HeatCapacity = 1;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 0;

f = fopen([cd,'\materials\hb-curves\','steel_1008.tab'],'r');
HB = fscanf(f,'%f');
fclose(f);
m.MagneticPermeability = [HB(1:2:end),HB(2:2:end)];

m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'steel_1008',m)

%% ===>  m15_26g
m.ThermalConductivity = 30;
m.HeatCapacity = 1;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 0;

f = fopen([cd,'\materials\hb-curves\','m15_26g.tab'],'r');
HB = fscanf(f,'%f');
fclose(f);
m.MagneticPermeability = [HB(1:2:end),HB(2:2:end)];

m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'m15_26g',m)

%% ===>  m530_50a
m.ThermalConductivity = 30;
m.HeatCapacity = 1;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 0;

f = fopen([cd,'\materials\hb-curves\','m530_50a.tab'],'r');
HB = fscanf(f,'%f');
fclose(f);
m.MagneticPermeability = [HB(1:2:end),HB(2:2:end)];

m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'m530_50a',m)

%% ===> material name: NdFe30
m.ThermalConductivity = 25;
m.HeatCapacity = 1;
m.ElectricPermitivity = pcts.e0;
m.ElectricConductivity = 0;
m.MagneticPermeability = 1.0445730167132*pcts.mu0;

m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'NdFe30',m)

%% aluminium
m.ThermalConductivity = 209;
m.HeatCapacity = 2;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 38e6;
m.MagneticPermeability = pcts.mu0;
m.MassDensity = 7800;
m.YoungModulus = 1;
m.PoissonRatio = 0.3;
WriteMDataBinary(fdir,'aluminium',m)

%% win
m.ThermalConductivity = 0.83;
m.HeatCapacity = 2;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 1;
m.MagneticPermeability = 4000*pcts.mu0;
m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'win',m)

%% win
m.ThermalConductivity = 0.83;
m.HeatCapacity = 2;
m.ElectricPermitivity = [pcts.e0 pcts.e0];
m.ElectricConductivity = 1;
m.MagneticPermeability = 4000*pcts.mu0;
m.MassDensity = 7800;
m.YoungModulus = 10;
m.PoissonRatio = 0.2;
WriteMDataBinary(fdir,'Y30BH',m)

end

