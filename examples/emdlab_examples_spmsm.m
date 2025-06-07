% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));
% design variables
gv_ISD = 74;
gv_OSD = 125;
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
gv_embrace = 0.89;
% define geometry data base
g = emdlab_g2d_db;
% mesh density function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2-gv_dm-gv_g,gv_ISD/2,gv_OSD/2], [3,0.5,0.5,3], r, 'linear','extrap');
% add geometry from library
emdlab_g2d_lib_stc1(g, gv_ISD, gv_OSD, gv_Ns, gv_wst, gv_dss, gv_bs0, gv_hs0, gv_tta);
emdlab_g2d_glib_spm1(g, gv_Dsh, gv_ISD-2*gv_g, gv_p, gv_dm, gv_embrace)
% set mesh density function
g.setMeshLengthByRadialFunction(f_mesh);
g.showSketch;
m = g.generateMesh('mg0');
% set mesh zone colors
m.setMeshZoneColor('stator',200,200,200)
m.setMeshZoneColor('sc',255,137,39);
m.setMeshZoneColor('rotor',200,200,200)
m.setMeshZoneColor('magnet',28,255,28)
m.setMeshZoneColor('sap',0,255,255)
m.setMeshZoneColor('rap',0,255,255)
% add materials
m.addMaterial('air');
m.addMaterial('es_M530_50A');
m.addMaterial('copper');
% set materials
m.setMaterial('rap','air');
m.setMaterial('sap','air');
m.setMaterial('rotor','es_M530_50A');
m.setMaterial('stator','es_M530_50A');
m.setMaterial('sc','copper');
m.setMaterial('magnet','air');
% generate full mesh
m.aux_cmrjmzx('stator',gv_Ns,gv_Ns)
m.aux_cmrjmzx('sap',gv_Ns,gv_Ns)
m.aux_cmrjmzx('rotor',gv_p,gv_p)
m.aux_cmrjmzx('rap',gv_p,gv_p)
m.aux_cmrjmzx('sc',gv_Ns,1)
m.aux_crmz('sc', gv_Ns, gv_Ns)
m.aux_cmrjmzx('magnet',gv_p,1)
m.aux_crmz('magnet', gv_p, gv_p)
% generate air gap mesh
m.ggmesh;
kr = m.getfbnioc([0,0],gv_ISD/2-gv_g);
ks = m.getfbnioc([0,0],gv_ISD/2);
rps = m.nodes(kr,:);
sps = m.nodes(ks,:);
% getting a moving contact object
agm = emdlab_mcs_circularAirGap(rps, sps, 4);
% adding air gap to mesh zones
m.addmz('AG', agm.m);
m.setMaterial('AG','air');
m.setMeshZoneColor('AG',153,217,234);
% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
% setting physical units
s.setUnit('length', 'mm');
s.setUnit('currentDensity', 'A/mm^2');
s.setUnit('magneticVectorPotential', 'A/m');
s.setDepth(70);
% set magnetization
Hc = -922100;
for i = 1:2:gv_p
s.setMagnetization(['magnet',num2str(i)],emdlab_solvers_ms2d_magnetization(Hc,'-r'));
end
for i = 2:2:gv_p
s.setMagnetization(['magnet',num2str(i)],emdlab_solvers_ms2d_magnetization(Hc,'r'));
end
% apply boundary conditions
m.ggmesh;
m.showmzs;
s.setAzBC(m.getfbn, 0);
% solve and plot results
figure;
s.setSolverRelativeError(1e-7);
s.solve
f = s.plotBmagF(18,'rotor', 'stator', 'magnet1', 'magnet2', 'magnet3', 'magnet4');
s.plotBrBtOnCircle(0,0,gv_ISD/2-gv_g/2, 1000);


