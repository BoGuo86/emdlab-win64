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

writeParFile('Rso',p.Rso,'Ts',p.Ts,'bss1',p.bss1,'Tsso1',p.Tsso1,...
    'wtb',p.wtb,'hss1',p.hss1,'alphass1',p.alphass1,'rsh',p.rsh,'wsy',p.wsy,...
    'rss1',p.rss1,'rss2',p.rss2,'g',0.5,'Rro',100,'mesh_thetag',1.7,...
    'mesh_thetas',7,'mesh_thetar',3);

m = TMDBC;
glib_outer_rotor_pmsm;
m.read_g2d_bin('geom.g2d','MM1');


m.cmmz('s2','s1',[1,0]);
for i = 1:2:2*(p.Ns-1)
    m.crmz(['s',num2str(i+2)],['s',num2str(i)],p.Ts);
    m.crmz(['s',num2str(i+3)],['s',num2str(i+1)],p.Ts);
end
temp = getlist('s',1:2*p.Ns);
m.jmzs('stator',temp{:});

m.setMaterial('c11','copper');
m.cmmz('c21','c11',[1,0]);

for i = 1:p.Ns-1
    m.crmz(['c1',num2str(i+1)],['c1',num2str(i)],p.Ts);
    m.crmz(['c2',num2str(i+1)],['c2',num2str(i)],p.Ts);
    m.crmz(['sap',num2str(i+1)],['sap',num2str(i)],p.Ts);
end

m.mzs.stator.moveNodes


theta = linspace(-pi/26,pi/26,14);
mp = zeros(length(theta),2);
for i = 1:length(theta)
    mp(i,1) = 92*cos(theta(i));
    mp(i,2) = 92*sin(theta(i));
end
zvec = linspace(0,7,7);
Power = 1.2;
zvec = zvec.^Power/7^(Power-1);
[a,b] = tmzpc_swipRadial0(mp,diff(zvec));
m.addmz('magnet1',TMZPC(b,a));

zvec = linspace(0,7,5);
Power = 1.2;
zvec = zvec.^Power/7^(Power-1);
[a,b] = tmzpc_swipRadial0(a(end-14+1:end,:),diff(zvec));
m.addmz('r1',TMZPC(b,a));

m.setmzc('magnet1','c');
for i = 1:p.Nm-1
    m.crmz(['magnet',num2str(i+1)],['magnet',num2str(i)],p.Tr);
    if rem(i,2)
        m.setmzc(['magnet',num2str(i)],'g');
    end
    m.crmz(['r',num2str(i+1)],['r',num2str(i)],p.Tr);
end
temp = getlist('r',1:p.Nm);
m.jmzs('rotor',temp{:});

m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'m19_24ga');
m.addMaterial(materialdir,'copper');
% setting materials
m.setMaterial('stator','m19_24ga');

m.setMaterial('rotor','m19_24ga');

ttt= 25.714299999999980;
m.rmz('stator',ttt);
for j = 1:p.Ns
    m.rmz(['c1',num2str(j)],ttt);
    m.rmz(['c2',num2str(j)],ttt);
    m.rmz(['sap',num2str(j)],ttt);
end
    
m.ggmesh;

kr = m.getnIndexOnCircle([0,0],p.Rso);
ks = m.getnIndexOnCircle([0,0],p.Rso+p.g);
rps = m.nodes(kr,:);
sps = m.nodes(ks,:);

agm = mc_circularAirGap(rps, sps);
%% Calling Solver
s = IHNLNRMSTL3(m);
s.setUnit('length', 'mm');
s.setUnit('currentDensity', 'A/mm^2');
s.setUnit('magneticVectorPotential', 'A/m');
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

s.defineMatrix('PhaseA', p.Ncoil, 1);
s.defineMatrix('PhaseB', p.Ncoil, 1);
s.defineMatrix('PhaseC', p.Ncoil, 1);

s.addMeshZone2Matrix('PhaseA', 'c11', 1);
s.addMeshZone2Matrix('PhaseA', 'c21', -1);
s.addMeshZone2Matrix('PhaseA', 'c12', -1);
s.addMeshZone2Matrix('PhaseA', 'c22', 1);
s.addMeshZone2Matrix('PhaseA', 'c16', 1);
s.addMeshZone2Matrix('PhaseA', 'c26', -1);
s.addMeshZone2Matrix('PhaseA', 'c110', -1);
s.addMeshZone2Matrix('PhaseA', 'c210', 1);
s.addMeshZone2Matrix('PhaseA', 'c114', 1);
s.addMeshZone2Matrix('PhaseA', 'c214', -1);
s.addMeshZone2Matrix('PhaseA', 'c118', -1);
s.addMeshZone2Matrix('PhaseA', 'c218', 1);
s.addMeshZone2Matrix('PhaseA', 'c119', 1);
s.addMeshZone2Matrix('PhaseA', 'c219', -1);

s.addMeshZone2Matrix('PhaseB', 'c15', 1);
s.addMeshZone2Matrix('PhaseB', 'c25', -1);
s.addMeshZone2Matrix('PhaseB', 'c14', -1);
s.addMeshZone2Matrix('PhaseB', 'c24', 1);
s.addMeshZone2Matrix('PhaseB', 'c18', 1);
s.addMeshZone2Matrix('PhaseB', 'c28', -1);
s.addMeshZone2Matrix('PhaseB', 'c19', -1);
s.addMeshZone2Matrix('PhaseB', 'c29', 1);
s.addMeshZone2Matrix('PhaseB', 'c113', 1);
s.addMeshZone2Matrix('PhaseB', 'c213', -1);
s.addMeshZone2Matrix('PhaseB', 'c117', -1);
s.addMeshZone2Matrix('PhaseB', 'c217', 1);
s.addMeshZone2Matrix('PhaseB', 'c121', 1);
s.addMeshZone2Matrix('PhaseB', 'c221', -1);

s.addMeshZone2Matrix('PhaseC', 'c17', 1);
s.addMeshZone2Matrix('PhaseC', 'c27', -1);
s.addMeshZone2Matrix('PhaseC', 'c13', -1);
s.addMeshZone2Matrix('PhaseC', 'c23', 1);
s.addMeshZone2Matrix('PhaseC', 'c112', 1);
s.addMeshZone2Matrix('PhaseC', 'c212', -1);
s.addMeshZone2Matrix('PhaseC', 'c111', -1);
s.addMeshZone2Matrix('PhaseC', 'c211', 1);
s.addMeshZone2Matrix('PhaseC', 'c115', 1);
s.addMeshZone2Matrix('PhaseC', 'c215', -1);
s.addMeshZone2Matrix('PhaseC', 'c116', -1);
s.addMeshZone2Matrix('PhaseC', 'c216', 1);
s.addMeshZone2Matrix('PhaseC', 'c120', 1);
s.addMeshZone2Matrix('PhaseC', 'c220', -1);


theta_old = 0;
%% proccess
rotorAngle = linspace(0,2*p.Tr,30);
Ntheta = length(rotorAngle);
LinkageFlux = zeros(3,Ntheta);
Bxs = zeros(s.m.mzs.stator.Ne,Ntheta);
Bys = zeros(s.m.mzs.stator.Ne,Ntheta);
Bxr = zeros(s.m.mzs.rotor.Ne,Ntheta);
Byr = zeros(s.m.mzs.rotor.Ne,Ntheta);

shifttpm = [0.1,0];

s.m.shmz('rotor',shifttpm);
    for j = 1:p.Nm
        s.m.shmz(['magnet',num2str(j)],shifttpm);
    end
    agm.shiftOuter(shifttpm);
    
% loop for calculation of linkage flux 
for i = 1:Ntheta
    % rotation of regions
    s.rotateMeshZone('rotor',rotorAngle(i)-theta_old, shifttpm);
    for j = 1:p.Nm
        s.rotateMeshZone(['magnet',num2str(j)],rotorAngle(i)-theta_old, shifttpm);
    end
    agm.rotateOuter(rotorAngle(i)-theta_old, shifttpm);
    s.addmz('AG',agm.m);
    theta_old = rotorAngle(i);
    s.m.ggmesh;
    % getting index fo boundary conditions
    k0 = s.m.edges(s.m.bedges,1:2);
    k0 = unique(k0(:));
    s.bcs.clearAllBCs;
    s.bcs.setDirichlet(k0,0);
    % calling and running solver
 
    s.solve;
    LinkageFlux(1,i) = s.evalMatrixLinkageFlux('PhaseA');
LinkageFlux(2,i) = s.evalMatrixLinkageFlux('PhaseB');
LinkageFlux(3,i) = s.evalMatrixLinkageFlux('PhaseC');

    Bxs(:,i) = s.results.Bex(s.m.ezi(:,s.m.mzs.stator.zi));
    Bys(:,i) = s.results.Bey(s.m.ezi(:,s.m.mzs.stator.zi));
    Bxr(:,i) = s.results.Bex(s.m.ezi(:,s.m.mzs.rotor.zi));
    Byr(:,i) = s.results.Bey(s.m.ezi(:,s.m.mzs.rotor.zi));
    s.removemz('AG');
end



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
