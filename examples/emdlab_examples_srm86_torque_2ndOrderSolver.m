% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));

Nposition = 30;
rotorPosition = linspace(0,30,Nposition);
Delta_Theta = (pi/6)/(Nposition-1);
torque_Arkkio = zeros(1,Nposition);
torque_MST = Nposition;
Wcoenergy = Nposition;

x = tic;
parfor i = 1:Nposition    
    [Wcoenergy(i),torque_MST(i)] = runMSA(rotorPosition(i));
end
clc;
toc(x);

% torque calculated previously by ANSYS Maxwell
fid = fopen('srm-te-maxwell.tab','r');
fgetl(fid);
data = fscanf(fid,'%f');
fclose(fid);
data = reshape(data,2,[])';
te = data(:,2);

% plot torque curves
figure;
title('Static Torque Curve');
hold on; box on;
plot(linspace(0,30,30),te,'marker','^');
% plot(rotorPosition,torque_Arkkio,'marker','s');
% plot(rotorPosition,torque_MST,'marker','o');
plot(rotorPosition,diff_cen(Wcoenergy,Delta_Theta),'marker','*');

function [Wcoenergy, torque_MST] = runMSA(rotorPosition)
% variables
gv_Ns = 8;
gv_Nr = 6;
gv_Dsh = 20;
gv_ISD = 74;
gv_OSD = 125;
gv_gap = 0.5;
gv_beta_r = 19.8 * pi/180;
gv_beta_s = 17.1 * pi/180;
gv_wry = 7.53;
gv_wsy = 7.7;
% mesh density function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2,gv_OSD/2], [gv_wry/1,0.5,gv_wsy/1], r);
% generation of geometry
g = emdlab_g2d_db;
emdlab_g2d_lib_srm1(g, gv_Ns, gv_Nr, gv_Dsh, gv_ISD, gv_OSD, gv_gap, gv_beta_r, gv_beta_s, gv_wry, gv_wsy)
% mesh generation
g.setMeshLengthByRadialFunction(f_mesh);
m = g.generateMesh('mg0');
m.setPrintFlag(0);
% set mesh zone colors
m.setMeshZoneColor('Rotor',90,90,90);
m.setMeshZoneColor('Stator',90,90,90);
m.setMeshZoneColor('RotorAP',200,191,231);
m.setMeshZoneColor('sca',255,137,39);
% add materials
m.addMaterial('air');
m.addMaterial('es_M530_50A');
m.addMaterial('copper');
% set mesh zone materials
m.setMaterial('RotorAP','air');
m.setMaterial('Rotor','es_M530_50A');
m.setMaterial('Stator','es_M530_50A');
m.setMaterial('sca','copper');
% construct full mesh
m.aux_cmrjmzx('Rotor', gv_Nr, gv_Nr);
m.aux_cmrjmzx('Stator', gv_Ns, gv_Ns);
m.aux_cmrjmzx('RotorAP', gv_Nr, gv_Nr);
m.aux_cmrmzx('sca', gv_Ns, gv_Ns);
% rotate moving mesh zones
m.rotateMeshZone('Rotor', (rotorPosition-30)*pi/180);
m.rotateMeshZone('RotorAP', (rotorPosition-30)*pi/180);
% generate air gap mesh
m.ggmesh;
kr = m.getfbnioc([0,0],gv_ISD/2-gv_gap);
ks = m.getfbnioc([0,0],gv_ISD/2);
rps = m.nodes(kr,:);
sps = m.nodes(ks,:);
% getting a moving contact object
agm = emdlab_mcs_circularAirGap(rps, sps, 1);
% adding air gap to mesh zones
m.addmz('AG', agm.m);
m.setMaterial('AG','air');
m.setMeshZoneColor('AG',153,217,234);
% getting an instance of solver object
m.ggmesh;
m.gd2elements;
s = emdlab_solvers_ms2d_tl6_ihnl(m);
% setting physical units
s.setUnit('length', 'mm');
s.setUnit('current', 'A');
s.setDepth(70);
% define new winding
s.defineWinding('winding_1');
s.addMeshZone2Winding('winding_1','sca11',278,'positive');
s.addMeshZone2Winding('winding_1','sca21',278,'negative');
s.addMeshZone2Winding('winding_1','sca15',278,'negative');
s.addMeshZone2Winding('winding_1','sca25',278,'positive');
% set winding current
s.setWindingCurrent('winding_1', 3.6);
% apply boundary conditions
m.ggmesh;
s.setAzBC(m.getfbn, 0);
% solve and plot results
s.setSolverRelativeError(1e-4)
s.setMonitor(0);
s.solve
s.setSolverMaxIteration(100);

[~,Wcoenergy] = s.evalTotalEnergyCoenergy;
torque_MST = s.evalTorqueByMST(gv_ISD/2-gv_gap/2);

end

