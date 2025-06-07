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

wnames = cell(1,nr*nd);
for i = 1:nd
    for j = 1:nr
        wnames{nr*(i-1)+j} = ['c',num2str(i),num2str(j)];
        g = g.newcface(wnames{nr*(i-1)+j},...
            [-Wslot/2+j*Wmid+(2*j-1)*Rwire,i*Wmid+(2*i-1)*Rwire],Rwire,'Nnodes',10);
    end
end

g = g.newdswdkps([-Wslot/2,0],[Wslot/2,0],'Nnodes',5*nr);
g = g.newdswdkps([Wslot/2,0],[Wslot/2,Hslot],'Nnodes',5*nd);
g = g.newdswdkps([Wslot/2,Hslot],[-Wslot/2,Hslot],'Nnodes',5*nr);
g = g.newdswdkps([-Wslot/2,Hslot],[-Wslot/2,0],'Nnodes',5*nd);
g = g.newloop('slot','S1',1,'S2',1,'S3',1,'S4',1);

g = g.newface('slot','slot',wnames{:});


m = TMDBC();
m = m.addmz('slot',g.getdmmz('slot'));
m = m.addMaterial(materialdir,'air');
for i = 1:nd
    for j = 1:nr
        m = m.addmz(wnames{nr*(i-1)+j},g.getdmmz(wnames{nr*(i-1)+j}));
    end
end
% m = m.ggmesh;

% %% Calling Solver
% s = IHNLNRMSTL3(m);clear m
% s.units.length = 'mm';
% s.units.currentDensity = 'A/mm^2';
% s.units.magneticVectorPotential = 'A/m';
% %% Process
% for i = 1:nd
%     for j = 1:nr
%         s = s.setExcitation(wnames{nr*(i-1)+j},1,'C');
%     end
% end
% k0 = s.m.getnIndexOnLine([-Wslot/2,Hslot],[Wslot/2,Hslot]);
% s = s.setdbc(k0,0);
% % calling and running solver
% s.m = s.m.evalKeFeC('TL3');
% s = s.assignEdata;
% s = s.solve(1e-2,20);
% s.plotAmag

