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


    
m.ggmesh;

kr = m.getnIndexOnCircle([0,0],p.Rso);
ks = m.getnIndexOnCircle([0,0],p.Rso+p.g);
rps = m.nodes(kr,:);
sps = m.nodes(ks,:);

agm = mc_circularAirGap(rps, sps);
m.addmz('AG',agm.m);
%% Calling Solver
s = IHNLNRMDWTMTL3(m);
s.setUnit('length', 'mm');
s.setUnit('currentDensity', 'A/mm^2');
s.setUnit('magneticVectorPotential', 'A/m');
s.setUnit('time', 'ms');
s.setDepth(p.Lst);
% %% Setting Magnets Excitations
% Hc = 3.5014e+05;
% for i = 1:p.Nm
%     if rem(i,2) == 0
%         s.setMagnetization(['magnet',num2str(i)],...
%             getRotate(msMagnetization(-Hc,'r'),i*p.Tr));
%     else
%         s.setMagnetization(['magnet',num2str(i)],...
%             getRotate(msMagnetization(Hc,'r'),i*p.Tr));
%     end
% end

s.defineWinding('PhaseA', 'wtype', 'stranded', 'extype', 'voltage','exvalue', tf_sin('amplitude',20,'frequency',50));
% s.defineWinding('PhaseB', p.Ncoil, 1);
% s.defineWinding('PhaseC', p.Ncoil, 1);

s.defineCoil('c11', 'turns', p.Ncoil, 'direction', 'positive');
s.addCoil2Winding('PhaseA', 'c11');
s.defineCoil('c21', 'turns', p.Ncoil, 'direction', 'negative');
s.addCoil2Winding('PhaseA', 'c21');
s.defineCoil('c12', 'turns', p.Ncoil, 'direction', 'negative');
s.addCoil2Winding('PhaseA', 'c12');
s.defineCoil('c22', 'turns', p.Ncoil, 'direction', 'positive');
s.addCoil2Winding('PhaseA', 'c22');
s.defineCoil('c16', 'turns', p.Ncoil, 'direction', 'positive');
s.addCoil2Winding('PhaseA', 'c16');
s.defineCoil('c26', 'turns', p.Ncoil, 'direction', 'negative');
s.addCoil2Winding('PhaseA', 'c26');
s.defineCoil('c110', 'turns', p.Ncoil, 'direction', 'negative');
s.addCoil2Winding('PhaseA', 'c110');
s.defineCoil('c210', 'turns', p.Ncoil, 'direction', 'positive');
s.addCoil2Winding('PhaseA', 'c210');
s.defineCoil('c114', 'turns', p.Ncoil, 'direction', 'positive');
s.addCoil2Winding('PhaseA', 'c114');
s.defineCoil('c214', 'turns', p.Ncoil, 'direction', 'negative');
s.addCoil2Winding('PhaseA', 'c214');
s.defineCoil('c118', 'turns', p.Ncoil, 'direction', 'negative');
s.addCoil2Winding('PhaseA', 'c118');
s.defineCoil('c218', 'turns', p.Ncoil, 'direction', 'positive');
s.addCoil2Winding('PhaseA', 'c218');
s.defineCoil('c119', 'turns', p.Ncoil, 'direction', 'positive');
s.addCoil2Winding('PhaseA', 'c119');
s.defineCoil('c219', 'turns', p.Ncoil, 'direction', 'negative');
s.addCoil2Winding('PhaseA', 'c219');

% s.addMeshZone2Matrix('PhaseB', 'c15', 1);
% s.addMeshZone2Matrix('PhaseB', 'c25', -1);
% s.addMeshZone2Matrix('PhaseB', 'c14', -1);
% s.addMeshZone2Matrix('PhaseB', 'c24', 1);
% s.addMeshZone2Matrix('PhaseB', 'c18', 1);
% s.addMeshZone2Matrix('PhaseB', 'c28', -1);
% s.addMeshZone2Matrix('PhaseB', 'c19', -1);
% s.addMeshZone2Matrix('PhaseB', 'c29', 1);
% s.addMeshZone2Matrix('PhaseB', 'c113', 1);
% s.addMeshZone2Matrix('PhaseB', 'c213', -1);
% s.addMeshZone2Matrix('PhaseB', 'c117', -1);
% s.addMeshZone2Matrix('PhaseB', 'c217', 1);
% s.addMeshZone2Matrix('PhaseB', 'c121', 1);
% s.addMeshZone2Matrix('PhaseB', 'c221', -1);
% 
% s.addMeshZone2Matrix('PhaseC', 'c17', 1);
% s.addMeshZone2Matrix('PhaseC', 'c27', -1);
% s.addMeshZone2Matrix('PhaseC', 'c13', -1);
% s.addMeshZone2Matrix('PhaseC', 'c23', 1);
% s.addMeshZone2Matrix('PhaseC', 'c112', 1);
% s.addMeshZone2Matrix('PhaseC', 'c212', -1);
% s.addMeshZone2Matrix('PhaseC', 'c111', -1);
% s.addMeshZone2Matrix('PhaseC', 'c211', 1);
% s.addMeshZone2Matrix('PhaseC', 'c115', 1);
% s.addMeshZone2Matrix('PhaseC', 'c215', -1);
% s.addMeshZone2Matrix('PhaseC', 'c116', -1);
% s.addMeshZone2Matrix('PhaseC', 'c216', 1);
% s.addMeshZone2Matrix('PhaseC', 'c120', 1);
% s.addMeshZone2Matrix('PhaseC', 'c220', -1);


%% Boundary condition

s.m.ggmesh;
s.bcs.clearAllBCs;
s.bcs.setDirichlet(s.m.getfbn,0);

s.setSimulationTime(1, 200);

s.saveBe(1)
s.solve
