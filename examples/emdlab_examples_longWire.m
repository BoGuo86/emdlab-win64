% the purpose of this example is to compare the first- and second-order solvers
% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));
% define geometry data base
g = emdlab_g2d_db;
% add loops
l1t = g.addCircleLoop(0,0,5);
l2t = g.addCircleLoop(0,0,50);
% define faces
g.addFace('wire', l1t);
g.addFace('air', l2t, l1t);
% set maximum mesh lngth
f_mesh = @(r) interp1([0,5,51],[1,3,15],r);
g.setMeshLengthByRadialFunction(f_mesh);
% g.showSketch;
m = g.generateMesh('mg0');
% set mesh zonecolor
m.setMeshZoneColor('wire',255,127,39);
m.setMeshZoneColor('air',150,237,239);
% add materials
m.addMaterial('air');
m.addMaterial('copper');
% set material of the regions
m.setMaterial('wire','air');
m.setMaterial('air','air');
% getting an instance of solver object
s = emdlab_solvers_ms2d_tl3_ihnl(m);
% setting physical units
s.setUnit('length', 'mm');
s.setUnit('current', 'A');
s.setDepth(1000);
% set boundary conditions
m.showmzs;
m.ggmesh;
s.defineWinding('wire');
s.addMeshZone2Winding('wire', 'wire', 1, "positive");
s.setWindingCurrent('wire',10000);
s.setAzBC(m.getfbn,0);
% run solver
s.setMonitor(0);
s.solve;
s.plotAmag
% plot By(x)
[bx, by] = s.getBvecOnSegment([-50,0],[50,0],1000);
figure;
plot(bx(:,1), by(:,2),'r', 'LineWidth',1);
% plot B mag
s.plotBmag
clim([0,0.4]);
% using second-order solver
m.gd2elements;
s2 = emdlab_solvers_ms2d_tl6_ihnl(m);
% setting physical units
s2.setUnit('length', 'mm');
s2.setUnit('current', 'A');
s2.setDepth(1000);
% set boundary conditions
s2.defineWinding('wire');
s2.addMeshZone2Winding('wire', 'wire', 1, "positive");
s2.setWindingCurrent('wire',10000);
s2.setAzBC(m.getfbn,0);
s2.setMonitor(0);
s2.solve;
s2.plotAmag;
s2.evalBe;
s2.plotBmag;
clim([0,0.4]);
figure(3);
x = linspace(-50,50,1000);
[bx,by] = s2.getBxByOnPoints(x,0*x);
hold on
plot(x,by,'b', 'LineWidth',1);
% analytic solution
figure(3);
x = linspace(-50,50,1000);
y = 0*x;
index = abs(x)<=5;
y(index) = 4*pi*1e-7*10000*x(index)*1e3/(2*pi*5^2);
y(~index) = 4*pi*1e-7*10000*1e3./(2*pi*x(~index));
plot(x,y,'k', 'LineWidth',1);
box on;
xlabel('$x$-axis [mm]', 'FontSize',14, 'Interpreter','latex');
ylabel('$B_y(x)$ [tesla]', 'FontSize',14, 'Interpreter','latex');
legend('First-order triangle', 'Second-order triangle', 'Exact solution');