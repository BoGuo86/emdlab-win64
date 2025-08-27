% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));
% variables
gv_Ns = 6;
gv_Nr = 4;
gv_Dsh = 20;
gv_ISD = 60;
gv_OSD = 125;
gv_gap = 0.45;
gv_beta_r = 34.2 * pi/180;
gv_beta_s = 30 * pi/180;
gv_wry = 10.42;
gv_wsy = 10;
rotorPosition = -20;
% dependent variables
gv_wrt = 2 * (gv_ISD/2-gv_gap) * sin(gv_beta_r/2);
gv_gammar = asin(gv_wrt*0.5/(gv_Dsh/2+gv_wry));
gv_alpharp = 2*pi/gv_Nr;
gv_wst = 2 * (gv_ISD/2) * sin(gv_beta_s/2);
gv_gammas = asin(gv_wst*0.5/(gv_OSD/2-gv_wsy));
gv_alphasp = 2*pi/gv_Ns;
% mesh density function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2,gv_OSD/2], [gv_wry/3,0.4,gv_wsy/3], r, 'linear', 'extrap');
% generation of geometry
g = emdlab_g2d_db;
% adding points
p1t = g.addPoint(0,0);
p2t = g.addPoint(gv_Dsh/2,0);
p3t = g.addPoint(gv_ISD/2 - gv_gap,0);
p4t = g.addPoint((gv_ISD/2 - gv_gap)*cos(gv_beta_r/2),(gv_ISD/2 - gv_gap)*sin(gv_beta_r/2));
p5t = g.addPoint((gv_Dsh/2+gv_wry)*cos(gv_gammar),(gv_Dsh/2+gv_wry)*sin(gv_gammar));
p6t = g.addPoint((gv_Dsh/2+gv_wry)*cos(gv_alpharp/2),(gv_Dsh/2+gv_wry)*sin(gv_alpharp/2));
p7t = g.addPoint((gv_Dsh/2)*cos(gv_alpharp/2),(gv_Dsh/2)*sin(gv_alpharp/2));
p8t = g.addPoint((gv_ISD/2 - gv_gap)*cos(gv_alpharp/2),(gv_ISD/2 - gv_gap)*sin(gv_alpharp/2));
p9t = g.addPoint(gv_ISD/2,0);
p10t = g.addPoint(gv_OSD/2,0);
p11t = g.addPoint((gv_OSD/2)*cos(gv_alphasp/2),(gv_OSD/2)*sin(gv_alphasp/2));
p12t = g.addPoint((gv_OSD/2-gv_wsy)*cos(gv_alphasp/2),(gv_OSD/2-gv_wsy)*sin(gv_alphasp/2));
p13t = g.addPoint((gv_OSD/2-gv_wsy)*cos(gv_gammas),(gv_OSD/2-gv_wsy)*sin(gv_gammas));
p14t = g.addPoint((gv_ISD/2)*cos(gv_beta_s/2),(gv_ISD/2)*sin(gv_beta_s/2));
p15t = g.addPoint((gv_ISD/2)*cos(gv_alphasp/2),(gv_ISD/2)*sin(gv_alphasp/2));
% adding edges
e1t = g.addSegment(p2t, p3t);
e2t = g.addArc(p1t, p3t, p4t, 1);
e3t = g.addSegment(p4t, p5t);
e4t = g.addArc(p1t, p5t, p6t, 1);
e5t = g.addSegment(p6t, p7t);
e6t = g.addArc(p1t, p7t, p2t, 0);
e7t = g.addArc(p1t, p4t, p8t, 1);
e8t = g.addSegment(p8t, p6t);
e9t = g.addSegment(p9t, p10t);
e10t = g.addArc(p1t, p10t, p11t, 1);
e11t = g.addSegment(p11t, p12t);
e12t = g.addArc(p1t, p12t, p13t, 0);
e13t = g.addSegment(p13t, p14t);
e14t = g.addArc(p1t, p14t, p9t, 0);
e15t = g.addArc(p1t, p15t', p14t, 0);
e16t = g.addSegment(p12t, p15t);
% adding loops
l1t = g.addLoop(e1t,e2t,e3t,e4t,e5t,e6t);
l2t = g.addLoop(-e3t,e7t,e8t,-e4t);
l3t = g.addLoop(e9t,e10t,e11t,e12t,e13t,e14t);
l4t = g.addLoop(-e13t,-e12t,e16t,e15t);
% adding faces
g.addFace('Stator', l3t);
g.addFace('Rotor', l1t);
g.addFace('RotorAP', l2t);
g.addFace('sca', l4t);
% mesh generation
g.setMeshLengthByRadialFunction(f_mesh);
m = g.generateMesh('mg0');
m.setPrintFlag(0);
% set mesh zone colors
m.setMeshZoneColor('Rotor',90,90,90);
m.setMeshZoneColor('Stator',90,90,90);
m.setMeshZoneColor('RotorAP',0,255,255);
m.setMeshZoneColor('sca',255,137,39);
% add materials
m.addMaterial('m530', emdlab_mlib_es_M530_50A);
m.addMaterial('copper', emdlab_mlib_copper);
% set materials
m.setMaterial('Rotor','m530');
m.setMaterial('Stator','m530');
m.setMaterial('sca','copper');
% construct full mesh
m.aux_cmxjcrj('Rotor', gv_Nr);
m.aux_cmxjcrj('Stator', gv_Ns);
m.aux_cmxjcrj('RotorAP', gv_Nr);
m.aux_cmxcr('sca', gv_Ns);
% rotate moving mesh zones
m.rotateMeshZone('Rotor', rotorPosition*pi/180);
m.rotateMeshZone('RotorAP', rotorPosition*pi/180);
% generate air gap mesh
m.aux_addCircularAirGap('ag',0,0,gv_ISD/2-gv_gap,0,0,gv_ISD/2,2)
% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
% set unit of the length
s.setLengthUnit('mm');
s.setDepth(70);
% define new winding
s.defineCoil('phaseA');
s.addMeshZone2Coil('phaseA','sca11',375,1);
s.addMeshZone2Coil('phaseA','sca21',375,-1);
s.addMeshZone2Coil('phaseA','sca14',375,-1);
s.addMeshZone2Coil('phaseA','sca24',375,1);
% set winding current
s.setCoilCurrent('phaseA', 3.6);
% apply boundary conditions
s.setAzBC(m.getfbn, 0);
% solve and plot results
s.setSolverRelativeError(1e-4);
s.setSolverMaxIteration(100);
s.solve;
s.gui;