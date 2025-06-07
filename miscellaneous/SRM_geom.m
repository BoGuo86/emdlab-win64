%% initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
p = srm0_design;
p.set_ns(12);
p.set_nr(8);
Ncoil = 400;
%% Creation of mesh and geometry
m = TMDBC;
p.writeParFile;
geom_srm0;
m.read_g2d_bin('geom.g2d');

m.setmzColor('s1',[133, 193, 233]/255);
m.setmzColor('r1',[133, 193, 233]/255);

m.setMaterial('r1','m19_24ga');
m.setMaterial('s1','m19_24ga');
m.setMaterial('c11','copper');
m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'m19_24ga');
m.addMaterial(materialdir,'copper');

m.cmirrormz('s2','s1',[1,0]);
for i = 1:2:2*(p.ns/p.Nsplit-1)
    m.crotatemz(['s',num2str(i+2)],['s',num2str(i)],p.Ts);
    m.crotatemz(['s',num2str(i+3)],['s',num2str(i+1)],p.Ts);
end
temp = getlist('s',1:2*p.ns/p.Nsplit);
m.joinmzs('stator',temp{:});

m.cmirrormz('c21','c11',[1,0]);
m.setmzColor('c11',[241,195,15]/255);
m.setmzColor('c21','c');
for i = 1:p.ns/p.Nsplit-1
    m.crotatemz(['c1',num2str(i+1)],['c1',num2str(i)],p.Ts);
    m.crotatemz(['c2',num2str(i+1)],['c2',num2str(i)],p.Ts);
end

m.cmirrormz('r2','r1',[cos(p.Tr/2),sin(p.Tr/2)]);
for i = 1:2:2*(p.nr/p.Nsplit-1)
    m.crotatemz(['r',num2str(i+2)],['r',num2str(i)],p.Tr);
    m.crotatemz(['r',num2str(i+3)],['r',num2str(i+1)],p.Tr);
end
temp = getlist('r',1:2*p.nr/p.Nsplit);
m.joinmzs('rotor',temp{:});

m.setmzColor('rap1',[195, 155, 211]/255)
for i = 1:p.nr/p.Nsplit-1
    m.crotatemz(['rap',num2str(i+1)],['rap',num2str(i)],p.Tr);
end

m.rotatemz('rotor',-p.Tr/2);
for i = 1:p.nr/p.Nsplit
    m.rotatemz(['rap',num2str(i)],-p.Tr/2);
end

m.ggmesh;
kr = m.getnIndexOnCircle([0,0],p.Rro);
ks = m.getnIndexOnCircle([0,0],p.Rro+p.g);
rps = m.nodes(kr,:);
sps = m.nodes(ks,:);

s = IHNLNRMSTL3(m);clear m
s.units.length = 'mm';
s.units.currentDensity = 'A/mm^2';
s.units.magneticVectorPotential = 'A/m';

pos_a = {'c11','c24','c17','c210'};
neg_a = {'c21','c14','c27','c110'};

%% proccess
phaseCurrent = linspace(0,20,20);
rotorAngle = linspace(0,p.Tr/2,2);
Ni = length(phaseCurrent);
Ntheta = length(rotorAngle);
linkageFlux = zeros(Ntheta,Ni);
theta_old = 0;
% loop for calculation of linkage flux 
for i = 1:Ntheta
    % rotation of regions
    s.m.rotatemz('rotor',rotorAngle(i)-theta_old);
        for j = 1:p.nr/p.Nsplit
            s.m.rotatemz(['rap',num2str(j)],rotorAngle(i)-theta_old);
        end
    if p.Nsplit>1
        s.addmz('AG',getmz(rotaterps(MC_AAG(sps,rps),rotorAngle(i))));
    else
        s.addmz('AG',getmz(rotaterps(MC_AG(sps,rps),rotorAngle(i))));
    end
    theta_old = rotorAngle(i);
    s.m.ggmesh;
    % getting index fo boundary conditions
    k0 = [s.m.getnIndexOnCircle([0,0],p.rsh);...
        s.m.getnIndexOnCircle([0,0],p.Rso)];
    s.clearallbcs;
    s.setdbc(k0,0);
    if p.Nsplit>1
        k = s.m.getfb;
        k = setdiff(unique(k(:)),k0);
        [km,ks] = s.m.splitPeriodic(k,2*pi/p.Nsplit);
        s.setopbc(km',ks');
    end
    % calling solver
    s.m.evalKeFe('TL3');
    for j = 1:Ni
        % set excitation
        for k = 1:4
            s.setExcitation(pos_a{k},phaseCurrent(j)*Ncoil/p.np,'C');
            s.setExcitation(neg_a{k},-phaseCurrent(j)*Ncoil/p.np,'C');
        end
        s.assignEdata;
        s.solve(1e-8,30);
        % eval linkage flux
        for k = 1:4
            linkageFlux(i,j) = linkageFlux(i,j) + s.evalLF(pos_a{k}) - s.evalLF(neg_a{k});
        end
    end
    s.removemz('AG');
end
% linkageFlux = linkageFlux*p.nser*p.Lst*Ncoil/p.np/1000;
close all
hold all
tmpC = [phaseCurrent,fliplr(phaseCurrent)];
tmpL = [linkageFlux(1,:),fliplr(linkageFlux(end,:))];
for i=1:Ntheta
    for j = 1:Ni
        plot(phaseCurrent,linkageFlux(i,:),...
            'color','b','marker','s','markersize',8,'linestyle','--',...
            'linewidth',1.5);
    end
end
% title(['Wcon = ',num2str(polyarea(tmpC,tmpL)),' Wcon,des = ',num2str(p.Wcon)]);
xlabel('Phase Current [A]');
ylabel('Linkage Flux [wb]');