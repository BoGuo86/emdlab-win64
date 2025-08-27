% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));
% design variables
gv_ISD = 75;
gv_OSD = 125;
gv_Lstk = 70;
gv_Ns = 36;
gv_p = 4;
gv_wst = 3;
gv_dss = 15;
gv_bs0 = 2.4;
gv_hs0 = 0.6;
gv_tta = 25;
gv_Dsh = 24;
gv_Ntc = 85;
gv_g = 0.4;
gv_Iph = 1.9;
gv_Nfb = 4;
% define geometry data base
g = emdlab_g2d_db;
% mesh density function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2,gv_OSD/2], [2,0.4,3], r, 'linear','extrap');
% add geometry from library
emdlab_g2d_lib_tc3(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta, 'stator', 'sc', 'sap');
emdlab_g2d_lib_rm_synrm1(g, gv_Dsh, gv_ISD-2*gv_g, gv_p, 0.45, ones(1,gv_Nfb)*0.8,  ones(1,gv_Nfb-1)*0.8, 0.6, 'rotor', 'fb')
% set mesh density function
g.setMeshLengthByRadialFunction(f_mesh);
% g.showSketch;
m = g.generateMesh('mg0');
if gv_Nfb>1
    m.joinMeshZones('fb','fb'+string(1:gv_Nfb));
else
    m.changeMeshZoneName('fb1', 'fb');
end
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
m.aux_cmxjcrj('fb',gv_p)
m.aux_cmxjcr('sc',gv_Ns)
% rotate mesh zones
for mz = ["rotor", "fb"]
    m.rotateMeshZone(mz,pi/4-pi/gv_Ns);
end
% generate air gap mesh
m.aux_addCircularAirGap('ag',0,0,gv_ISD/2-gv_g,0,0,gv_ISD/2,2);
% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
s.setLengthUnit('mm');
s.setDepth(gv_Lstk);
% define windings
s.defineCoil('phaseA');
s.defineCoil('phaseB');
s.defineCoil('phaseC');
% index of each coil assocciating with each winding
pAn = [1,2,3]; pAn = [pAn, pAn+18];
pAp = [10,11,12]; pAp = [pAp, pAp+18];
pBp = pAp + 6;
pBn = pAn + 6;
pCp = pAn + 3;
pCn = pAp + 3;
% assingation of mesh zones to windings
for i = 1:6
    s.addMeshZone2Coil('phaseA', ['sc',num2str(pAp(i))], gv_Ntc, 1);
    s.addMeshZone2Coil('phaseA', ['sc',num2str(pAn(i))], gv_Ntc, -1);
    s.addMeshZone2Coil('phaseB', ['sc',num2str(pBp(i))], gv_Ntc, 1);
    s.addMeshZone2Coil('phaseB', ['sc',num2str(pBn(i))], gv_Ntc, -1);
    s.addMeshZone2Coil('phaseC', ['sc',num2str(pCp(i))], gv_Ntc, 1);
    s.addMeshZone2Coil('phaseC', ['sc',num2str(pCn(i))], gv_Ntc, -1);
end
% set phase currents
s.setCoilCurrent('phaseA', gv_Iph*1.41);
s.setCoilCurrent('phaseB', -gv_Iph*1.41/2);
s.setCoilCurrent('phaseC', -gv_Iph*1.41/2);
% apply boundary conditions
s.setAzBC(m.getfbn, 0);
% solve and plot results
s.setSolverRelativeError(1e-4);
s.setSolverMaxIteration(100);
s.solve;
s.gui;
