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
geom_outer_rotor_pmsm;
m.read_g2d_bin('geom.g2d');

m.setmzColor('s1',[133, 193, 233]/255);
m.cmirrormz('s2','s1',[1,0]);
for i = 1:2:2*(p.Ns-1)
    m.crotatemz(['s',num2str(i+2)],['s',num2str(i)],p.Ts);
    m.crotatemz(['s',num2str(i+3)],['s',num2str(i+1)],p.Ts);
end
temp = getlist('s',1:2*p.Ns);
m.joinmzs('stator',temp{:});

m.setMaterial('c11','copper');
m.cmirrormz('c21','c11',[1,0]);
m.setmzColor('c11',[241,195,15]/255);
m.setmzColor('c21','c');
m.setmzColor('sap1','w');
for i = 1:p.Ns-1
    m.crotatemz(['c1',num2str(i+1)],['c1',num2str(i)],p.Ts);
    m.crotatemz(['c2',num2str(i+1)],['c2',num2str(i)],p.Ts);
    m.crotatemz(['sap',num2str(i+1)],['sap',num2str(i)],p.Ts);
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

m.setmzColor('magnet1','c');
for i = 1:p.Nm-1
    m.crotatemz(['magnet',num2str(i+1)],['magnet',num2str(i)],p.Tr);
    if rem(i,2)
        m.setmzColor(['magnet',num2str(i)],'g');
    end
    m.crotatemz(['r',num2str(i+1)],['r',num2str(i)],p.Tr);
end
temp = getlist('r',1:p.Nm);
m.joinmzs('rotor',temp{:});
m.setmzColor('rotor',[133, 193, 233]/255);

m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'m19_24ga');
m.addMaterial(materialdir,'copper');
% setting materials
m.setMaterial('stator','m19_24ga');

m.setMaterial('rotor','m19_24ga');

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
%% Setting Magnets Excitations
% Hc = 3e5;
% for i = 1:p.Nm
%     if rem(i,2) == 0
%         s.setMagnetization(['magnet',num2str(i)],...
%             getRotate(msMagnetization(Hc,'r'),i*p.Tr));
%     else
%         s.setMagnetization(['magnet',num2str(i)],...
%             getRotate(msMagnetization(-Hc,'r'),i*p.Tr));
%     end
% end
% theta_old = 0;
%% proccess
LinkageFlux = zeros(3,3);
ia = [1,0,0];
ib = [0,1,0];
ic = [0,0,1];


pos_a = {'c11','c22','c16','c210','c114','c119','c218'};
neg_a = {'c21','c12','c26','c110','c118','c214','c219'};

pos_b = {'c15','c24','c18','c29','c113','c217','c121'};
neg_b = {'c25','c14','c28','c19','c117','c213','c221'};

pos_c = {'c17','c23','c112','c211','c115','c120','c216'};
neg_c = {'c212','c27','c13','c111','c116','c215','c220'};

s.addmz('AG',getmz(MC_AG(sps,rps)));

s.m.ggmesh;

% getting index fo boundary conditions
    k0 = s.m.edges(s.m.bedges,1:2);
    k0 = unique(k0(:));
    s.clearallbcs;
    s.setdbc(k0,0);
    % calling and running solver
    s.m.evalKeFe('TL3');

% loop for calculation of linkage flux 
for i = 1:3
    
    
    
    
    for j  =1:7
        s.setExcitation(pos_a{j},p.Ncoil*ia(i),'C');
        s.setExcitation(neg_a{j},-p.Ncoil*ia(i),'C');
        s.setExcitation(pos_b{j},p.Ncoil*ib(i),'C');
        s.setExcitation(neg_b{j},-p.Ncoil*ib(i),'C');
        s.setExcitation(pos_c{j},p.Ncoil*ic(i),'C');
        s.setExcitation(neg_c{j},-p.Ncoil*ic(i),'C');
    end
    
    s.assignEdata;
    s.solve(1e-2,20);

for j  =1:7
    LinkageFlux(1,i) = LinkageFlux(1,i) + s.evalLF(pos_a{j}) - s.evalLF(neg_a{j});
    LinkageFlux(2,i) = LinkageFlux(2,i) + s.evalLF(pos_b{j}) - s.evalLF(neg_b{j});
    LinkageFlux(3,i) = LinkageFlux(3,i) + s.evalLF(pos_c{j}) - s.evalLF(neg_c{j});
end


end

LinkageFlux = LinkageFlux*p.Ncoil*p.Lst/1000




