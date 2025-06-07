%% Initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
%% Inputs
% number of row conductors
nr = 1;
% number of height conductors
nd = 9;
% wire diameter
Dwire = 2.5;
% mid distance between wires
Wmid = 0.4;
% % wire radius
% gdata.Rwire = Dwire/2;
% minimum width of slot
Wslot = nr*(Dwire+2*Wmid);
% minimum height of slot
Hslot = nd*(Dwire+2*Wmid);
% stack length of slot
Lst = 400;
%% Creation of geometry and mesh
m = TMDBC;
geom_RectangularSlot(Dwire,Wmid);
m.read_g2d_bin('geom.g2d');
m.setmzColor('wire',[255, 128, 0]/255);
m.setmzColor('ins',[133, 193, 233]/255);
m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'copper');
m.setMaterial('wire','copper');
tmp = Wmid+Dwire/2;
m.cmirrormz('ins1','ins',[tmp,0],[tmp,1]);
m.cmirrormz('ins2','ins',[0,tmp],[1,tmp]);
m.cmirrormz('ins3','ins1',[0,tmp],[1,tmp]);
m.joinmzs('ins11','ins','ins1','ins2','ins3');
m.cmirrormz('wire1','wire',[tmp,0],[tmp,1]);
m.cmirrormz('wire2','wire',[0,tmp],[1,tmp]);
m.cmirrormz('wire3','wire1',[0,tmp],[1,tmp]);
m.joinmzs('c11','wire','wire1','wire2','wire3');
m.mzs.c11.moveNodes;
tmp = 2*Wmid+Dwire;
for i = 2:nd
    m.cshiftmz(['c',num2str(i),'1'],['c',num2str(i-1),'1'],[0,tmp]);
    m.cshiftmz(['ins',num2str(i),'1'],['ins',num2str(i-1),'1'],[0,tmp]);
end
for i = 1:nd
    for j = 2:nr
        m.cshiftmz(['c',num2str(i),num2str(j)],...
            ['c',num2str(i),num2str(j-1)],[tmp,0]);
         m.cshiftmz(['ins',num2str(i),num2str(j)],...
            ['ins',num2str(i),num2str(j-1)],[tmp,0]);
    end
end
enames = cell(1,nr*nd);
cnames = cell(1,nr*nd);
for i = 1:nd
    for j = 1:nr
        enames{nr*(i-1)+j} = ['ins',num2str(i),num2str(j)];
        cnames{nr*(i-1)+j} = ['c',num2str(i),num2str(j)];
    end
end
m.joinmzs('slot',enames{:});
m.ggmesh;
%% Calling Solver
s = IHLECTL3(m); clear m;
s.depth = Lst;
s.units.length = 'mm';
s.units.currentDensity = 'A/mm^2';
s.units.magneticVectorPotential = 'A/m';
s.units.frequency = 'Hz';
%% Process
s.defWinding('win','solid','current','series',sqrt(2),0,1);
for i = 1:nr
    
end
for i = 1:nd
    for j = 1:nr
        s.defCoil(['c',num2str(i),num2str(j)]);
        s.addCoil2Winding('win',['c',num2str(i),num2str(j)]);
    end
end
k0 = s.m.getnIndexOnLine([0,Hslot],[1,Hslot]);
s.clearallbcs;
s.setdbc(k0,0,0);
% calling and running solver
s.m.evalKeFe('TL3');
s.m.evalMe();
s.assignEdata;
%% Calling solver
frange = 0:10:100;
loss = zeros(1,length(frange));
for i = 1:length(frange)
    s.frequency = frange(i);
    s.solve;
    loss(i) = s.evalSolidLoss;
end
DCLoss = nr*nd*Lst/58e3/s.m.mzs.c11.area;
plot(frange,loss/DCLoss,'Marker','o')
title('Resistance Factor, kR = R_{AC}/R_{DC}');
xlabel('frequency [Hz]')