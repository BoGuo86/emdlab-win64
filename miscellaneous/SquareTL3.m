clc
clear
materialdir = [cd,'\MaterialsData'];
g = GDBC2D;

g = g.newpolygonalcb('s',[0,0;1,0;1,1;0,1]);
g = g.newd('s',inf,'s');

m = MDBCT(g);clear g;
m = m.addMaterial(materialdir,'air');
s = IHLTHSTL3(m);
s.m = s.m.ggmesh;

Nrefine = 6;
EE = zeros(1,Nrefine-1);

for i = 1:Nrefine
    
    if i >1
        s.m = s.m.strefine;
    end
    
    k1 = s.m.getIndexOnRay([0,0],[1,0]);
    k2 = s.m.getIndexOnRay([1,0],[1,1]);
    k3 = s.m.getIndexOnRay([1,1],[0,1]);
    k4 = s.m.getIndexOnRay([0,0],[0,1]);

    bc = @(p) sin(pi*p(:,2));
    s = s.setdbc(k2,0);
    bc = @(p) sin(pi*p(:,1)); 
    s = s.setdbc(k3,0);
    bc = @(p) sin(pi*p(:,2));
    s = s.setdbc(k4,0);
    bc = @(p) sin(pi*p(:,1)); 
    s = s.setdbc(k1,0);
    s = s.setExcitation('s',10);
    
    s.m = s.m.evalKeFeC('TL3');
    s = s.assignEdata;
    s = s.solve;
    s = s.evalQ;

    s = s.evalEtot;
    E(i) = s.Etot;
    if i>1
        s.Q = s.Q - [Q;Q;Q;Q];
        s = s.evalEtot;
        EE(i-1) = sqrt(s.Etot);
    end
    
    s = s.evalQ;
    Q = s.Q;
    
    
    
end

clc
EE(1:end-1)./EE(2:end)

s.plotT

