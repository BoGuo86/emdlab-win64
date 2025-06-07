% Initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
%% Mesh
m = TMDBC;
m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'copper');
m.addMaterial(materialdir,'m19_24ga');
m.addMaterial(materialdir,'iron');
% Rotor
Nr = 6;
Ns = 36;

Nsplit = gcd(Ns,Nr);

Nrsplit = Nr/Nsplit;
Nssplit = Ns/Nsplit;

glib_ipm_rt2(25,35,Nr,25,2,0.3,0.3,4,1,1.15);
m.read_g2d_bin('geom.g2d','MG1');
m.setMaterial('r1','m19_24ga');
m.setmzColor('r1',[0.7883,0.8529,0.4856]);
m.cmirrormz('r2','r1',[1,0]);
for i = 1:2:2*(Nrsplit-1)
    m.crotatemz(['r',num2str(i+2)],['r',num2str(i)],2*pi/Nr);
    m.crotatemz(['r',num2str(i+3)],['r',num2str(i+1)],2*pi/Nr);
end
temp = getlist('r',1:2*Nrsplit);
m.joinmzs('rotor',temp{:});
m.setmzColor('rap11',[0.9700,0.6951,0.2733]);

m.cmirrormz('rap21','rap11',[1,0]);
for i = 1:Nrsplit-1
    m.crotatemz(['rap1',num2str(i+1)],['rap1',num2str(i)],2*pi/Nr);
    m.crotatemz(['rap2',num2str(i+1)],['rap2',num2str(i)],2*pi/Nr);
end
m.setmzColor('m1',[0.3320,0.7487,0.6444]);
m.cmirrormz('m2','m1',[1,0]);
m.joinmzs('magnet1','m1','m2');
for i = 1:Nrsplit-1
    m.crotatemz(['magnet',num2str(i+1)],['magnet',num2str(i)],2*pi/Nr);
end

% stator

glib_im_st3(35.5,50,Ns,3,4,2,1,1,3.5,1.4);
m.read_g2d_bin('geom.g2d','MG1');
m.setMaterial('s1','m19_24ga');
m.setmzColor('s1',[0.3453,0.9468,0.5202]);

m.cmirrormz('s2','s1',[1,0]);

for i = 1:2:2*(Nssplit-1)
    m.crotatemz(['s',num2str(i+2)],['s',num2str(i)],2*pi/Ns);
    m.crotatemz(['s',num2str(i+3)],['s',num2str(i+1)],2*pi/Ns);
end
temp = getlist('s',1:2*Nssplit);
m.joinmzs('stator',temp{:});
m.setmzColor('c11',[0.8260,0.0210,0.0464]);
m.cmirrormz('c21','c11',[1,0]);
m.setmzColor('c21',[0.8925,0.3807,0.9237]);
for i = 1:(Nssplit-1)
    m.crotatemz(['c1',num2str(i+1)],['c1',num2str(i)],2*pi/Ns);
    m.crotatemz(['c2',num2str(i+1)],['c2',num2str(i)],2*pi/Ns);
    m.crotatemz(['sap',num2str(i+1)],['sap',num2str(i)],2*pi/Ns);
end

m.ggmesh

kr = m.getnIndexOnCircle([0,0],35);
ks = m.getnIndexOnCircle([0,0],35.5);
rps = m.nodes(kr,:);
sps = m.nodes(ks,:);
m.addmz('AG',getmz(MC_AAG(sps,rps)));

%% Calling Solver
s = IHNLNRMSTL3(m);
s.units.length = 'mm';
s.units.currentDensity = 'A/mm^2';
s.units.magneticVectorPotential = 'A/m';

%% Setting Magnets Excitations
Hc = 1e6;
for i = 1:Nrsplit
    if rem(i,2) == 0
        s.setMagnetization(['magnet',num2str(i)],...
            getRotate(msMagnetization(Hc,'x'),(i-1)*2*pi/Nr));
    else
        s.setMagnetization(['magnet',num2str(i)],...
            getRotate(msMagnetization(-Hc,'x'),(i-1)*2*pi/Nr));
    end
end

s.m.ggmesh;
k0 = [s.m.getnIndexOnCircle([0,0],25);s.m.getnIndexOnCircle([0,0],50)];
k = s.m.getfbn;
k = setdiff(k,k0);
[km,ks] = s.m.splitPeriodic(k,2*pi/Nsplit);
s.clearallbcs;
s.setdbc(k0,0);
if rem(Nrsplit,2)==0
s.setepbc(km',ks');
else
    s.setopbc(km',ks');
end

% calling and running solver

s.solve(1e-6,20);
s.plotBrBtOnCircle([0,0],35.25)
