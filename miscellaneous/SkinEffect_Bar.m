%% Initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
%% Inputs
Wbar = 10;
Hbar = 20;
nw = 10;
nh = 20;
%% Creation of Mesh
[x,y] = meshgrid(linspace(-Wbar/2,Wbar/2,nw),linspace(0,Hbar,nh));
m = TMDBC;
m.addmz('bar',TMZPC(delaunay([x(:),y(:)]),[x(:),y(:)]));
m.addMaterial(materialdir,'copper');
m.setMaterial('bar','copper');
m.ggmesh;
%% Calling Solver
s = IHLECTL3(m);clear m;
s.setUnit('length', 'mm');
s.setUnit('currentDensity', 'A/mm^2');
s.setUnit('magneticVectorPotential', 'A/m');
s.setUnit('frequency', 'Hz');
s.setDepth(1000);
%% Process
s.defineCoil('bar');
s.defineWinding('bar','wtype','solid','extype','current','sptype','series','amplitude',1000,'phase',0);
s.addCoil2Winding('bar','bar');
k0 = s.m.getnIndexOnLine([-Wbar/2,Hbar],[Wbar/2,Hbar]);
s.bcs.clearAllBCs
s.bcs.setDirichlet(k0,0);
%% Calling solver
frange = 0:5:100;
loss = zeros(1,length(frange));
for i = 1:length(frange)
    s.setGlobalFrequency(frange(i));
    s.solve;
    loss(i) = s.evalSolidLoss;
end
DCLoss = s.getDepth/58e6/s.m.mzs.bar.getArea/1e-6;
plot(frange,loss/DCLoss,'Marker','o')
title('Resistance Factor, kR = R_{AC}/R_{DC}');
xlabel('frequency [Hz]')
