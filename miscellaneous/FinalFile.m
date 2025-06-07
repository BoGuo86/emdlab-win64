% Initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
%% Machine Data
p = pmsm0_design;
p.set_Pout(0.5);
p.set_Rso(70);
p.set_Ns(24);
p.set_Nm(8);
p.set_wtb(6);
p.set_wsy(5);
p.set_bss1(1);
p.set_hss1(1);
p.set_D(0.65);
p.set_betam(0.85);
p.set_thetas(2);
p.set_thetag(2);
p.set_thetar(4);
p.set_Span(3);
%% Creation of Mesh
m = TMDBC;
p.writeParFile;
geom_pmsm0;
m.read_g2d_bin('geom.g2d','MG0');
% adding needed materials
m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'m19_24ga');
m.addMaterial(materialdir,'copper');
% setting materials
m.setMaterial('s1','m19_24ga');
m.setMaterial('c11','copper');
m.setMaterial('r1','m19_24ga');
% set colors
m.setmzColor('rap1','w');
m.setmzColor('sap1','w');
m.setmzColor('s1',[133, 193, 233]/255);
m.setmzColor('r1',[133, 193, 233]/255);
m.cmirrormz('s2','s1',[1,0]);
for i = 1:2:2*(p.Ns-1)
    m.crotatemz(['s',num2str(i+2)],['s',num2str(i)],p.Ts);
    m.crotatemz(['s',num2str(i+3)],['s',num2str(i+1)],p.Ts);
end
m.cmirrormz('c21','c11',[1,0]);
m.setmzColor('c11',[241,195,15]/255);
m.setmzColor('c21','c');
for i = 1:p.Ns-1
    m.crotatemz(['c1',num2str(i+1)],['c1',num2str(i)],p.Ts);
    m.crotatemz(['c2',num2str(i+1)],['c2',num2str(i)],p.Ts);
    m.crotatemz(['sap',num2str(i+1)],['sap',num2str(i)],p.Ts);
end
m.cmirrormz('r2','r1',[1,0]);
for i = 1:2:2*(p.Nm-1)
    m.crotatemz(['r',num2str(i+2)],['r',num2str(i)],p.Tr);
    m.crotatemz(['r',num2str(i+3)],['r',num2str(i+1)],p.Tr);
end
m.cmirrormz('m2','m1',[1,0]);
m.joinmzs('magnet1','m1','m2');
m.setmzColor('magnet1',[195, 155, 211]/255)
for i = 1:p.Nm-1
    m.crotatemz(['magnet',num2str(i+1)],['magnet',num2str(i)],p.Tr);
    m.crotatemz(['rap',num2str(i+1)],['rap',num2str(i)],p.Tr);
end
temp = getlist('s',1:2*p.Ns);
m.joinmzs('stator',temp{:});
temp = getlist('r',1:2*p.Nm);
m.joinmzs('rotor',temp{:});
m.ggmesh;
kr = m.getnIndexOnCircle([0,0],p.Rro);
ks = m.getnIndexOnCircle([0,0],p.Rro+p.g);
rps = m.nodes(kr,:);
sps = m.nodes(ks,:);
%% Calling Solver
s = IHNLNRMSTL3(m);
s.units.length = 'mm';
s.units.currentDensity = 'A/mm^2';
s.units.magneticVectorPotential = 'A/m';
%% Setting Magnets Excitations
Hc = 1e6;
for i = 1:p.Nm
    if rem(i,2) == 0
        s.setMagnetization(['magnet',num2str(i)],...
            getRotate(msMagnetization(Hc,'r'),i*p.Tr));
    else
        s.setMagnetization(['magnet',num2str(i)],...
            getRotate(msMagnetization(-Hc,'r'),i*p.Tr));
    end
end
theta_old = 0;
%% proccess
rotorAngle = linspace(0,p.Tr,15);
Ntheta = length(rotorAngle);
Bxs = zeros(s.m.mzs.stator.Ne,Ntheta);
Bys = zeros(s.m.mzs.stator.Ne,Ntheta);
Bxr = zeros(s.m.mzs.rotor.Ne,Ntheta);
Byr = zeros(s.m.mzs.rotor.Ne,Ntheta);
% loop for calculation of linkage flux 
for i = 1:Ntheta
    % rotation of regions
    s.m.rotatemz('rotor',rotorAngle(i)-theta_old);
    for j = 1:p.Nm
        s.m.rotatemz(['rap',num2str(j)],rotorAngle(i)-theta_old);
        s.m.rotatemz(['magnet',num2str(j)],rotorAngle(i)-theta_old);
    end
    s.addmz('AG',getmz(rotaterps(MC_AG(sps,rps),rotorAngle(i))));
    theta_old = rotorAngle(i);
    s.m.ggmesh;
    % getting index fo boundary conditions
    k0 = s.m.edges(s.m.bedges,1:2);
    k0 = unique(k0(:));
    s.clearallbcs;
    s.setdbc(k0,0);
    % calling and running solver
    s.m.evalKeFe('TL3');
    s.assignEdata;
    s.solve(1e-2,20);
    Bxs(:,i) = s.B(s.m.ezi(:,s.m.mzs.stator.zi),1);
    Bys(:,i) = s.B(s.m.ezi(:,s.m.mzs.stator.zi),2);
    Bxr(:,i) = s.B(s.m.ezi(:,s.m.mzs.rotor.zi),1);
    Byr(:,i) = s.B(s.m.ezi(:,s.m.mzs.rotor.zi),2);
    s.removemz('AG');
end