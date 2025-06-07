%% Initialization
clc
clear
materialdir = [cd,'\MaterialsData'];

p.Rso = 183/2;
p.Ns = 21;
p.Ts = 2*pi/p.Ns;
p.bss1 = 5;
p.Tsso1 = 2*asin(p.bss1/2/p.Rso);
p.wtb = 5.2;
p.hss1 = 3;
p.alphass1 = 70;
p.rsh = 86/2;
p.wsy = 13.5;
p.rss1 = 4;
p.rss2 = 2;
p.Nm = 26;
p.Tr = 2*pi/p.Nm;
p.g = 0.5;
p.rpm = 450;
p.wm = p.rpm*pi/30;

myVar = 2/3;
p.Lst = 32*myVar;
p.Ncoil = 140*1/myVar;

m = TMDBC;
% glib_outer_rotor_pmsm;
m.read_m2df_bin('m2d\pmsm_outer');
% add materials
m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'m19_24ga');
m.addMaterial(materialdir,'iron');
% setting materials
m.setMaterial('stator','m19_24ga');
m.setMaterial('rotor','iron');

ttt= 25.714299999999980;
m.rotatemz('stator',ttt);
for j = 1:p.Ns
    m.rotatemz(['c1',num2str(j)],ttt);
    m.rotatemz(['c2',num2str(j)],ttt);
    m.rotatemz(['sap',num2str(j)],ttt);
end
    
m.ggmesh;

kr = m.getnIndexOnCircle([0,0],p.Rso);
ks = m.getnIndexOnCircle([0,0],p.Rso+p.g);
rps = m.nodes(kr,:);
sps = m.nodes(ks,:);

%% Calling Solver
s = IHNLNRMSTL3(m);
s.units.length = 'mm';
s.units.currentDensity = 'A/mm^2';
s.units.magneticVectorPotential = 'A/m';
s.setDepth(p.Lst);
%% Setting Magnets Excitations
Hc = 3.5014e+05;
for i = 1:p.Nm
    if rem(i,2) == 0
        s.setMagnetization(['magnet',num2str(i)],...
            getRotate(msMagnetization(-Hc,'r'),i*p.Tr));
    else
        s.setMagnetization(['magnet',num2str(i)],...
            getRotate(msMagnetization(Hc,'r'),i*p.Tr));
    end
end
theta_old = 0;
%% proccess
rotorAngle = linspace(0,2*p.Tr,30);
Ntheta = length(rotorAngle);
LinkageFlux = zeros(3,Ntheta);
Bxs = zeros(s.m.mzs.stator.Ne,Ntheta);
Bys = zeros(s.m.mzs.stator.Ne,Ntheta);
Bxr = zeros(s.m.mzs.rotor.Ne,Ntheta);
Byr = zeros(s.m.mzs.rotor.Ne,Ntheta);
% loop for calculation of linkage flux 
for i = 1:Ntheta
    % rotation of regions
    s.m.rotatemz('rotor',rotorAngle(i)-theta_old);
    for j = 1:p.Nm
        s.m.rotatemz(['magnet',num2str(j)],rotorAngle(i)-theta_old);
    end
    s.addmz('AG',getmz(rotatesps(MC_AG(sps,rps),rotorAngle(i))));
    theta_old = rotorAngle(i);
    s.m.ggmesh;
    % getting index fo boundary conditions
    k0 = s.m.edges(s.m.bedges,1:2);
    k0 = unique(k0(:));
    s.clearallbcs;
    s.setdbc(k0,0);
    % calling and running solver
 
    s.solve(1e-2,20);
%     s.plotBmag;
%     F(i) = getframe(gcf);
%     close(gcf);
LinkageFlux(1,i) = s.evalLF('c11')-s.evalLF('c21')-s.evalLF('c12')+s.evalLF('c22')...
    +s.evalLF('c16')-s.evalLF('c26')-s.evalLF('c110')+s.evalLF('c210')...
    +s.evalLF('c114')-s.evalLF('c214')-s.evalLF('c118')+s.evalLF('c218')...
    +s.evalLF('c119')-s.evalLF('c219');
LinkageFlux(2,i) = s.evalLF('c15')-s.evalLF('c25')-s.evalLF('c14')+s.evalLF('c24')...
    +s.evalLF('c18')-s.evalLF('c28')-s.evalLF('c19')+s.evalLF('c29')...
    +s.evalLF('c113')-s.evalLF('c213')-s.evalLF('c117')+s.evalLF('c217')...
    +s.evalLF('c121')-s.evalLF('c221');
LinkageFlux(3,i) = s.evalLF('c17')-s.evalLF('c27')-s.evalLF('c13')+s.evalLF('c23')...
    +s.evalLF('c112')-s.evalLF('c212')-s.evalLF('c111')+s.evalLF('c211')...
    +s.evalLF('c115')-s.evalLF('c215')-s.evalLF('c116')+s.evalLF('c216')...
    +s.evalLF('c120')-s.evalLF('c220');

    Bxs(:,i) = s.B(s.m.ezi(:,s.m.mzs.stator.zi),1);
    Bys(:,i) = s.B(s.m.ezi(:,s.m.mzs.stator.zi),2);
    Bxr(:,i) = s.B(s.m.ezi(:,s.m.mzs.rotor.zi),1);
    Byr(:,i) = s.B(s.m.ezi(:,s.m.mzs.rotor.zi),2);
    s.removemz('AG');
end

LinkageFlux = LinkageFlux*p.Ncoil;

BEMFs = p.wm*diff_cen(LinkageFlux,rotorAngle(2));

subplot(211)
hold all
plot(rotorAngle*13*180/pi,BEMFs(1,:),'marker','*');
plot(rotorAngle*13*180/pi,BEMFs(2,:),'marker','o');
plot(rotorAngle*13*180/pi,BEMFs(3,:),'marker','d');
set(gca,'xlim',[0,360])
xlabel('Electrical Degree')
ylabel('Volt')
legend({'E_a','E_b','E_c'})
title('Phase BEMFs')

subplot(212)
hold all
plot(rotorAngle*13*180/pi,BEMFs(1,:)-BEMFs(2,:),'marker','*');
plot(rotorAngle*13*180/pi,BEMFs(2,:)-BEMFs(3,:),'marker','o');
plot(rotorAngle*13*180/pi,BEMFs(3,:)-BEMFs(1,:),'marker','d');
set(gca,'xlim',[0,360])
xlabel('Electrical Degree')
ylabel('Volt')
legend({'E_{ab}','E_{b }','E_{ca}'})
title('Line to Line BEMFs')
