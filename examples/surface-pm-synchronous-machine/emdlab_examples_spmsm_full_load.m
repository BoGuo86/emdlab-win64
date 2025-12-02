% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));
% design variables
gv_ISD = 74;
gv_OSD = 125;
gv_Lstk = 70;
gv_Ns = 36;
gv_p = 4;
gv_wst = 3;
gv_dss = 15;
gv_bs0 = 2.4;
gv_hs0 = 0.6;
gv_tta = 25;
gv_Dsh = 28;
gv_dm = 3;
gv_g = 1;
gv_Ntc = 65;
gv_Is = 2.3;
gv_embrace = 0.82;
gv_Hc = -922100;
% define geometry data base
g = emdlab_g2d_db;
% mesh density function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2-gv_dm-gv_g,gv_ISD/2,gv_OSD/2], [4,1,0.5,3], r, 'linear','extrap');
% add geometry from library
emdlab_g2d_lib_tc3(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_rm_spm1(g, gv_Dsh, gv_ISD-2*gv_g, gv_p, gv_dm, gv_embrace, 'rotor', 'magnet', 'rap')
% set mesh density function
g.setMeshLengthByRadialFunction(f_mesh);
m = g.generateMesh('mm');
% add materials
m.addMaterial('m330', emdlab_mlib_es_M330_35A);
m.addMaterial('copper', emdlab_mlib_copper);
% set materials
m.setMaterial('rotor','m330');
m.setMaterial('stator','m330');
m.setMaterial('sc','copper');
% generate full mesh
m.aux_cmxjcrj('stator',gv_Ns)
m.aux_cmxjcrj('sap',gv_Ns)
m.aux_cmxjcrj('rotor',gv_p)
m.aux_cmxjcrj('rap',gv_p)
m.aux_cmxjcr('sc',gv_Ns)
m.aux_cmxjcr('magnet',gv_p)
% generate air gap mesh
m.aux_addCircularAirGapInterface('ag',0,0,gv_ISD/2-gv_g,0,0,gv_ISD/2,2,'inner')
% getting an instance of solver object
s = emdlab_solvers_mt2d_tl3_ihnlwm(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);
% set magnetization
for i = 1:2:gv_p
    s.setMagnetization(['magnet',num2str(i)],gv_Hc,'r');
    m.setMeshZoneColor(['magnet',num2str(i)],0,125,223);
end
for i = 2:2:gv_p
    s.setMagnetization(['magnet',num2str(i)],gv_Hc,'-r');
    m.setMeshZoneColor(['magnet',num2str(i)],255,70,70);
end
% define windings
s.defineCoil('phaseA','current','stranded');
s.defineCoil('phaseB','current','stranded');
s.defineCoil('phaseC','current','stranded');
% index of each coil assocciating with each winding
pAn = [1,2,3]; pAn = [pAn, pAn+18];
pAp = [10,11,12]; pAp = [pAp, pAp+18];
pBp = pAp + 6;
pBn = pAn + 6;
pCp = pAn + 3;
pCn = pAp + 3;
% assingation of mesh zones to windings
for i = 1:6
    s.addMeshZone2Coil('phaseA', ['sc',num2str(pAp(i))], gv_Ntc, 1, 0.15);
    s.addMeshZone2Coil('phaseA', ['sc',num2str(pAn(i))], gv_Ntc, -1, 0.15);
    s.addMeshZone2Coil('phaseB', ['sc',num2str(pBp(i))], gv_Ntc, 1, 0.15);
    s.addMeshZone2Coil('phaseB', ['sc',num2str(pBn(i))], gv_Ntc, -1, 0.15);
    s.addMeshZone2Coil('phaseC', ['sc',num2str(pCp(i))], gv_Ntc, 1, 0.15);
    s.addMeshZone2Coil('phaseC', ['sc',num2str(pCn(i))], gv_Ntc, -1, 0.15);
end
% define moving region
s.defineMovingRegion('moving_1', ["rotor", "rap", "magnet"+string(1:gv_p)], 'ag');
% apply boundary conditions
s.setAzBC(m.getfbn, 0);
% disable solver monitor
s.setMonitor(0);
% synchronous speed
wm = 1500*pi/30;
% phase currents in stator reference frame
alpha_Is_rrf = 90 * pi/180;
alpha_d4A = -110 * pi/180;
alpha_Is_srf = @(t) alpha_Is_rrf + alpha_d4A + 2*wm*t;
id_srf = @(t) gv_Is * cos(alpha_Is_srf(t));
iq_srf = @(t) gv_Is * sin(alpha_Is_srf(t));
ia_srf = @(t) id_srf(t);
ib_srf = @(t) -id_srf(t)/2 + sqrt(3)*iq_srf(t)/2;
ic_srf = @(t) -id_srf(t)/2 - sqrt(3)*iq_srf(t)/2;
s.setCoilCurrent('phaseA', ia_srf);
s.setCoilCurrent('phaseB', ib_srf);
s.setCoilCurrent('phaseC', ic_srf);
s.setSolverRelativeError(1e-3);
% total time steps for 20e-3 simulation
Nt = 200; 
% simulation time step size
timeStep = 20e-3/Nt;
% allocate memory to store calculated electric torque
te = zeros(1,Nt);
te_mst = zeros(1,Nt);
s.solveForInitialConditions;
for i = 1:Nt
    % calculate torque
    te(i) = s.evalTorqueByArkkio('ag', gv_g);
    te_mst(i) = s.evalTorqueByMST3(0,0,gv_ISD/2-gv_g/2,gv_g,1000);
    % rotate moving region
    s.rotateMovingRegion('moving_1', wm*timeStep);
    s.solveForOneTimeStep(timeStep);
end
% plot results
s.plotCoilCurrents;
s.plotCoilFluxLinkages;
s.plotCoilInducedVoltages;

id_srf = (s.coils.phaseA.current-0.5*s.coils.phaseB.current-0.5*s.coils.phaseC.current)*2/3;
iq_srf = (s.coils.phaseB.current-s.coils.phaseC.current)/sqrt(3);

fld_srf = (s.coils.phaseA.fluxLinkage-0.5*s.coils.phaseB.fluxLinkage-0.5*s.coils.phaseC.fluxLinkage)*2/3;
flq_srf = (s.coils.phaseB.fluxLinkage-s.coils.phaseC.fluxLinkage)/sqrt(3);

Te_dq = (3/2) * (4/2) * (fld_srf.*iq_srf - flq_srf.*id_srf);

figure;
hold on; box on;
plot(s.simTime, [te(1),te]);
plot(s.simTime, [te_mst(1),te_mst]);
plot(s.simTime, Te_dq);
xlabel('Time [s]'); ylabel('Electric torque [Nm]');
legend('Arkkio', 'MST', 'Te_{dq}');
set(gca,'ylim',[0,8])
hold off;


