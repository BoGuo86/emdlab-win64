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
gv_Iph = 1.9;
% define geometry data base
g = emdlab_g2d_db;
% mesh density function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2-gv_g,gv_OSD/2], [2,0.5,2], r, 'linear','extrap');
% add geometry from library
emdlab_g2d_lib_tc3(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_stta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_tc4(g, gv_Dsh, gv_ISD-2*gv_g, gv_Nr, gv_wrt, gv_drs, gv_br0, gv_hr0, gv_rtta, 'rotor', 'rbar', 'rap');
% set mesh density function
g.setMeshLengthByRadialFunction(f_mesh);
m = g.generateMesh('mg0');
m.joinMeshZones('rb','rbar','rap');
m.setMeshZoneColor('rb',100,100,100);
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
m.aux_cmxjcrj('stator',gv_Ns/2,2*pi/gv_Ns)
m.aux_cmxjcrj('sap',gv_Ns/2,2*pi/gv_Ns)
m.aux_cmxjcr('sc',gv_Ns/2,2*pi/gv_Ns)
m.aux_cmxjcrj('rotor',gv_Nr/2,2*pi/gv_Nr)
m.aux_cmxjcr('rb',gv_Nr/2,2*pi/gv_Nr)
% rotate rotor and bars
for mzName = ['rotor', 'rb' + string(1:13)]
    m.rotateMeshZone(mzName,pi/gv_Nr-pi/gv_Ns);
end
% add circular air gap
m.aux_addArcAirGap('ag',0,0,gv_ISD/2-gv_g,gv_ISD/2,2);
% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);
% define windings
s.defineCoil('phaseA');
s.defineCoil('phaseB');
s.defineCoil('phaseC');
% index of each coil assocciating with each winding
pAn = [1,2,3];
pAp = [10,11,12];
pBp = pAp + 6;
pBn = pAn + 6;
pCp = pAn + 3;
pCn = pAp + 3;
% assingation of mesh zones to windings
s.addMeshZone2Coil('phaseA', 'sc'+string(pAp), gv_Ntc, 1);
s.addMeshZone2Coil('phaseA', 'sc'+string(pAn), gv_Ntc, -1);
s.addMeshZone2Coil('phaseB', 'sc'+string(pBp), gv_Ntc, 1);
s.addMeshZone2Coil('phaseB', 'sc'+string(pBn), gv_Ntc, -1);
s.addMeshZone2Coil('phaseC', 'sc'+string(pCp), gv_Ntc, 1);
s.addMeshZone2Coil('phaseC', 'sc'+string(pCn), gv_Ntc, -1);
% set phase currents
s.setCoilCurrent('phaseA', gv_Iph*1.41);
s.setCoilCurrent('phaseB', -gv_Iph*1.41/2);
s.setCoilCurrent('phaseC', -gv_Iph*1.41/2);
% apply boundary conditions
[km,ks] = m.splitPeriodic(m.getfbniol_p0p1([0,0],[cos(-pi/gv_Ns),sin(-pi/gv_Ns)]), pi);
s.setAzBC(setdiff(m.getfbn,[km,ks]),0);
s.setEvenPeriodicBC(km,ks);
% solve and plot results
s.setSolverRelativeError(1e-5);
s.solve
s.gui;