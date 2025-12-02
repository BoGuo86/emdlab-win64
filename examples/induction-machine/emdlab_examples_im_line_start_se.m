% transient voltage-fed line-start simulation of induction motor
% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));
% design variables [mm]
gv_ISD = 74;
gv_OSD = 125;
gv_Lstk = 70;
gv_g = 0.5;
gv_Dsh = 24;
gv_Ns = 36;
gv_Nr = 26;
gv_wst = 3;
gv_dss = 15;
gv_bs0 = 2.4;
gv_hs0 = 0.5;
gv_stta = 35;
gv_wrt = 4;
gv_drs = 10;
gv_br0 = 2;
gv_hr0 = 0.7;
gv_rtta = 35;
gv_Ntc = 85;
gv_Vs = 325;
% define geometry data base
g = emdlab_g2d_db;
% mesh density function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2-gv_g,gv_OSD/2], [3,gv_g,3], r, 'linear','extrap');
% add geometry from library
emdlab_g2d_lib_tc3(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_stta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_tc4(g, gv_Dsh, gv_ISD-2*gv_g, gv_Nr, gv_wrt, gv_drs, gv_br0, gv_hr0, gv_rtta, 'rotor', 'rbar', 'rap');
% set mesh density function
g.setMeshLengthByRadialFunction(f_mesh);
m = g.generateMesh('mm');
m.joinMeshZones('rb','rbar','rap');
m.setMeshZoneColor('rb',0,192,232);
% add materials
m.addMaterial('m530', emdlab_mlib_es_M530_50A);
m.addMaterial('copper', emdlab_mlib_copper);
m.addMaterial('aluminium', emdlab_mlib_aluminium);
% set materials
m.setMaterial('rotor','m530');
m.setMaterial('stator','m530');
m.setMaterial('sc','copper');
m.setMaterial('rb','aluminium');
% generate full mesh
m.aux_cmxjcrj('stator',gv_Ns)
m.aux_cmxjcrj('sap',gv_Ns)
m.aux_cmxjcr('sc',gv_Ns)
m.aux_cmxjcrj('rotor',gv_Nr)
m.aux_cmxjcr('rb',gv_Nr)
m.shiftMeshZones(["rotor", "rb" + string(1:gv_Nr)], gv_g/2, 0);
% add circular air gap
m.aux_addCircularAirGapInterface('ag',gv_g/2,0,gv_ISD/2-gv_g,0,0,gv_ISD/2,1);
% getting an instance of solver object
s = emdlab_solvers_mt2d_tl3_ihnlwm(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);
% define windings
s.defineCoil('phaseA','voltage','stranded');
s.defineCoil('phaseB','voltage','stranded');
s.defineCoil('phaseC','voltage','stranded');
% index of each coil assocciating with each winding
pA = [1,2,3,10,11,12]; pA = [pA, pA+18];
ntc = gv_Ntc*[-1,-1,-1,1,1,1]; ntc = [ntc,ntc];
pB = pA + 6;
pC = pA + 3;
% assingation of mesh zones to windings
s.addMeshZones2Coil('phaseA', 'sc' + string(pA), ntc, 0.15);
s.addMeshZones2Coil('phaseB', 'sc' + string(pB), ntc, 0.15);
s.addMeshZones2Coil('phaseC', 'sc' + string(pC), -ntc, 0.15);
% set phase voltages
s.setCoilVoltage('phaseA', @(t) gv_Vs*sin(2*pi*50*t));
s.setCoilVoltage('phaseB', @(t) gv_Vs*sin(2*pi*50*t-2*pi/3));
s.setCoilVoltage('phaseC', @(t) gv_Vs*sin(2*pi*50*t+2*pi/3));
% define star connection
s.defineStarConnection('y1', 'phaseA', 'phaseB', 'phaseC');
% define rotor cage
s.defineCage('cage_1','rb'+string(1:gv_Nr));
% define moving region
s.defineMovingRegion('moving_1', ["rotor", "rb"+string(1:gv_Nr)], 'ag');
% apply boundary conditions
s.setAzBC(m.getfbn, 0);
% disable solver monitor
s.setMonitor(0);

% two vectors to store torque and speed
te = 0;
te_mst = 0;
wm = 0;

% simulation time step size
timeStep = 1e-3;

% run simulation for 400 time steps
for i = 1:400

    s.solveForOneTimeStep(timeStep);

    % calculate torque
%     te(end+1) = s.evalTorqueByArkkio('ag', gv_g);
    te_mst(end+1) = s.evalTorqueByMST3(gv_g/2,0,gv_ISD/2-gv_g/4,gv_g/2);

    % solve mechanical equation
    wm(end+1) = wm(end) + timeStep * (te_mst(end) - 5.1*sign(wm(end)) - 0.002*wm(end)) / 0.0025;

    % rotate moving region
    s.rotateMovingRegion('moving_1', wm(end)*timeStep, gv_g/2, 0);

end

s.plotCoilCurrents;

figure
plot(s.simTime, wm*30/pi)
xlabel('Time [s]'); ylabel('Rotor speed [rpm]');

id_srf = (s.coils.phaseA.current-0.5*s.coils.phaseB.current-0.5*s.coils.phaseC.current)*2/3;
iq_srf = (s.coils.phaseB.current-s.coils.phaseC.current)/sqrt(3);

fld_srf = (s.coils.phaseA.fluxLinkage-0.5*s.coils.phaseB.fluxLinkage-0.5*s.coils.phaseC.fluxLinkage)*2/3;
flq_srf = (s.coils.phaseB.fluxLinkage-s.coils.phaseC.fluxLinkage)/sqrt(3);

Te_dq = (3/2) * (4/2) * (fld_srf.*iq_srf - flq_srf.*id_srf);

figure;
hold on; box on;
plot(s.simTime, te);
plot(s.simTime, te_mst);
plot(s.simTime, Te_dq);
xlabel('Time [s]'); ylabel('Electric torque [Nm]');
legend('Arkkio', 'MST', 'Te_{dq}');
hold off;
