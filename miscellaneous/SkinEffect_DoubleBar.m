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
m.cshiftmz('bar2', 'bar', [0, Hbar]);
m.ggmesh;
%% Calling Solver
s = IHLECTL3(m);clear m;
s.units.length = 'mm';
s.units.currentDensity = 'A/mm^2';
s.units.magneticVectorPotential = 'A/m';
s.units.frequency = 'Hz';
s.setDepth(1000);
%% Process
s.defCoil('bar', 1);
s.defCoil('bar2', 1);
s.defWinding('bar','stranded','current','series',1000,0);
s.defWinding('bar2','stranded','current','series',1000,0);
s.addCoil2Winding('bar','bar');
s.addCoil2Winding('bar2','bar2');
k0 = s.m.getnIndexOnLine([-Wbar/2,2*Hbar],[Wbar/2,2*Hbar]);
s.clearallbcs;
s.setdbc(k0,0,0);
%% Calling solver
frange = 60;
loss = zeros(1,length(frange));
for i = 1:length(frange)
    s.setGlobalFrequency(frange(i));
    s.solve;
    loss(i) = s.evalSolidLoss
end
DCLoss = s.depth/58e6/s.m.mzs.bar.getArea/1e-3;
% plot(frange,loss/DCLoss,'Marker','o')
% title('Resistance Factor, kR = R_{AC}/R_{DC}');
% xlabel('frequency [Hz]')
s.evalStrandedLoss
