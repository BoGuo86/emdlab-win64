clc
clear
materialdir = [cd,'\MaterialsData'];
g = G2DKERNEL;

nr = 4;
nd = 7;
Rwire = 0.3;
Wmid = 0.12;
Wslot = 2*Rwire*nr+(nr+1)*Wmid;
Hslot = 2*Rwire*(nd)+2*nd*Wmid;

Center = [-Wslot/2+Wmid+Rwire,Wmid+Rwire];
g = g.newdarccppwdkps(Center,...
    [-Wslot/2+Wmid+2*Rwire,Wmid+Rwire],[-Wslot/2+Wmid+Rwire,Wmid+2*Rwire],'Nnodes',5);
g = g.newdarccppwdkps(Center,...
    [-Wslot/2+Wmid+Rwire,Wmid+2*Rwire],[-Wslot/2+Wmid,Wmid+Rwire],'Nnodes',5);
g = g.newdarccppwdkps(Center,...
    [-Wslot/2+Wmid,Wmid+Rwire],[-Wslot/2+Wmid+Rwire,Wmid],'Nnodes',5);
g = g.newdarccppwdkps(Center,...
    [-Wslot/2+Wmid+Rwire,Wmid],[-Wslot/2+Wmid+2*Rwire,Wmid+Rwire],'Nnodes',5);
g = g.newdswdkps([-Wslot/2,0],[-Wslot/2+Wmid+Rwire,0],'Nnodes',5);
g = g.newdswdkps([-Wslot/2+Wmid+Rwire,0],[-Wslot/2+Wmid+Rwire,Wmid],'Nnodes',4);
g = g.newdswdkps([-Wslot/2+Wmid,Wmid+Rwire],[-Wslot/2,Wmid+Rwire],'Nnodes',4);
g = g.newdswdkps([-Wslot/2,Wmid+Rwire],[-Wslot/2,0],'Nnodes',5);
g = g.newloop('c11','A1',1,'A2',1,'A3',1,'A4',1);
g = g.newloop('b1','S1',1,'S2',1,'A3',-1,'S3',1,'S4',1);
g = g.newface('c11','c11');
g = g.newface('b1','b1');

m = TMDBC();
m = m.addmz('c11',g.getdmmz('c11'));
m = m.addmz('b1',g.getdmmz('b1'));
m = m.cmirrormz('b2','b1',[-Wslot/2+Wmid+Rwire,0],[-Wslot/2+Wmid+Rwire,1]);
m = m.joinmzs('b3','b2','b1');
m = m.cmirrormz('b4','b3',[0,Wmid+Rwire],[1,Wmid+Rwire]);
m = m.joinmzs('e11','b4','b3');

for i = 2:nd
    m = m.cshiftmz(['c',num2str(i),'1'],['c',num2str(i-1),'1'],[0,2*Wmid+2*Rwire]);
    m = m.cshiftmz(['e',num2str(i),'1'],['e',num2str(i-1),'1'],[0,2*Wmid+2*Rwire]);
end
for i = 1:nd
    for j = 2:nr
        m = m.cshiftmz(['c',num2str(i),num2str(j)],...
            ['c',num2str(i),num2str(j-1)],[2*Wmid+2*Rwire,0]);
         m = m.cshiftmz(['e',num2str(i),num2str(j)],...
            ['e',num2str(i),num2str(j-1)],[2*Wmid+2*Rwire,0]);
    end
end
enames = cell(1,nr*nd);
cnames = cell(1,nr*nd);
for i = 1:nd
    for j = 1:nr
        enames{nr*(i-1)+j} = ['e',num2str(i),num2str(j)];
        cnames{nr*(i-1)+j} = ['c',num2str(i),num2str(j)];
    end
end
m = m.joinmzs('slot',enames{:});

for i = 1:nd
    for j = 1:nr
        m = m.setMaterial(cnames{nr*(i-1)+j},'copper');
    end
end

% m.showmzs
m = m.addMaterial(materialdir,'air');
m = m.addMaterial(materialdir,'copper');
m = m.ggmesh;

%% Calling Solver
s = IHNLNRMSTL3(m);clear m
s.units.length = 'mm';
s.units.currentDensity = 'A/mm^2';
s.units.magneticVectorPotential = 'A/m';
%% Process
for i = 1:nd
    for j = 1:nr
        s = s.setExcitation(['c',num2str(i),num2str(j)],1,'C');
    end
end
k0 = s.m.getnIndexOnLine([-Wslot/2,Hslot],[Wslot/2,Hslot]);
s = s.setdbc(k0,0);
% calling and running solver
s.m = s.m.evalKeFeC('TL3');
s = s.assignEdata;
s = s.solve(1e-2,20);

s.plotBmagSmooth

