clc
clear
materialdir = [cd,'\MaterialsData'];

g = GDBC2D;

p1 = getp('line',[0,0],[0.5,0],10,1.2);
p11 = getp('line',[0.5,0],[1,0],10,0.8);
p2 = getp('line',[1,0],[1,1],15,1.2);
p3 = getp('line',[1,1],[0,1],6,1);
p4 = [p2(:,1)-1,p2(:,2)];
p4 = preverse(p4);

p = [p1(1:end-1,:)
    p11(1:end-1,:)
    p2(1:end-1,:)
    p3(1:end-1,:)
    p4(1:end-1,:)];

g = g.newpolygonalcb('s',p);

g = g.newdDM('s','s');


m = MDBCT(g);
m = m.addMaterial(materialdir,'air');
s = IHLTHSTL3(m);
s.m = s.m.ggmesh;

Nrefine = 5;
EE = zeros(1,Nrefine-1);

for i = 1:Nrefine
    s.m = s.m.strefine;
    k1 = s.m.getIndexOnRay([0,0],[1,0]);
    k2 = s.m.getIndexOnRay([1,0],[1,1]);
    k3 = s.m.getIndexOnRay([1,1],[0,1]);
    k4 = s.m.getIndexOnRay([0,0],[0,1]);
    s = s.setdbc(k2,0);
    s = s.setdbc(k3,0);
    s = s.setdbc(k4,0);
    s = s.setdbc(k1,1);
%     s = s.setExcitation('s',10);
    s.m = s.m.evalKeFeC('TL3');
    s = s.assignEdata;
    s = s.solve;
    s = s.evalQ;
    
    if i>1
        s.Q = s.Q - [Q;Q;Q;Q];
        s = s.evalEtot;
        EE(i-1) = sqrt(s.Etot);
    end
    
    s = s.evalQ;
    Q = s.Q;
    
end

clc
for i = 2:length(EE)
    disp(EE(i-1)/EE(i))
end

s.plotT