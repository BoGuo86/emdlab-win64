% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));
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
% dependent variables
gv_wrt = 2 * (gv_ISD/2-gv_gap) * sin(gv_beta_r/2);
gv_gammar = asin(gv_wrt*0.5/(gv_Dsh/2+gv_wry));
gv_alpharp = 2*pi/gv_Nr;
gv_wst = 2 * (gv_ISD/2) * sin(gv_beta_s/2);
gv_gammas = asin(gv_wst*0.5/(gv_OSD/2-gv_wsy));
gv_alphasp = 2*pi/gv_Ns;
% mesh density function
f_mesh = @(r) interp1([gv_Dsh/2,gv_ISD/2,gv_OSD/2], [gv_wry/1,0.5,gv_wsy/1], r);
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
l1t = g.addLoop(e1t,1,e2t,1,e3t,1,e4t,1,e5t,1,e6t,1);
l2t = g.addLoop(e3t,0,e7t,1,e8t,1,e4t,0);
l3t = g.addLoop(e9t,1,e10t,1,e11t,1,e12t,1,e13t,1,e14t,1);
l4t = g.addLoop(e13t,0,e12t,0,e16t,1,e15t,1);
% adding faces
g.addFace('Rotor', l1t);
g.addFace('RotorAP', l2t);
g.addFace('Stator', l3t);
g.addFace('sca', l4t);
% mesh generation
% g.setMeshMaxLength(2);
g.setMeshLengthByRadialFunction(f_mesh);
g.showSketch;
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
m.rotateMeshZone('Rotor', -30*pi/180);
m.rotateMeshZone('RotorAP', -30*pi/180);
% generate air gap mesh
m.ggmesh;
kr = m.getfbnioc([0,0],gv_ISD/2-gv_gap);
ks = m.getfbnioc([0,0],gv_ISD/2);
rps = m.nodes(kr,:);
sps = m.nodes(ks,:);
% getting a moving contact object
agm = emdlab_mcs_circularAirGap(rps, sps, 2);
% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
% setting physical units
s.setUnit('length', 'mm');
s.setUnit('currentDensity', 'A/mm^2');
s.setUnit('magneticVectorPotential', 'A/m');
s.setDepth(70);
% define new winding
s.defineWinding('winding_1');
s.addMeshZone2Winding('winding_1','sca11',278,'positive');
s.addMeshZone2Winding('winding_1','sca21',278,'negative');
s.addMeshZone2Winding('winding_1','sca15',278,'negative');
s.addMeshZone2Winding('winding_1','sca25',278,'positive');
% set winding current
s.setWindingCurrent('winding_1', 3.6);
% number of rotor positions for calculation of the static torque curve
Nposition = 30;
rotorPosition = linspace(-30,0,Nposition);
Delta_Theta = (pi/6)/(Nposition-1);
torque_Arkkio = zeros(1,Nposition);
torque_MST = Nposition;
Wcoenergy = Nposition;
figure;
% loop for torque calculation
for i = 1:Nposition    

    % adding air gap to mesh zones    
    s.addmz('AG', agm.m);    
    s.m.setMaterial('AG','air');
    s.m.setMeshZoneColor('AG',153,217,234);

    % apply boundary conditions
    s.m.ggmesh;
    s.bcs.clearAllBCs;
    s.setAzBC(m.getfbn, 0);
    
    s.setMonitor(0);
    s.setSolverRelativeError(1e-7);
    s.solve
    s.setSolverMaxIteration(100);
    %s.plotBmagF(16,'Rotor', 'Stator');

    % torque calculation using different methods
    torque_Arkkio(i) = s.evalTorqueByArkkio('AG', gv_gap);
    torque_MST(i) = s.evalTorqueByMST([0,0], gv_ISD/2-gv_gap/2);
    [~,Wcoenergy(i)] = s.evalTotalEnergyCoenergy;

    % remove air-gap to remesh and add it again
    s.removemz('AG');

    % rotate moving mesh zones
    s.m.rotateMeshZone('Rotor', Delta_Theta);
    s.m.rotateMeshZone('RotorAP', Delta_Theta);
    agm.rotateInner(Delta_Theta);

end

% plot results
rotorPosition = rotorPosition + 30;

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
plot(rotorPosition,te,'marker','^');
plot(rotorPosition,torque_Arkkio,'marker','s');
plot(rotorPosition,torque_MST,'marker','o');
plot(rotorPosition,diff_cen(Wcoenergy,Delta_Theta),'marker','*');
xlabel('Rotor Position [deg]');
ylabel('Electric Torque [Nm]');
legend('ANSYS Maxwell', 'Arkkio', 'MST', 'Coenergy')

