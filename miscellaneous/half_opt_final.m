 alpha = 3.5
 
%% initialization
materialdir = [cd,'\MaterialsData'];
%% Mesh
m = TMDBC;
m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'m19_24ga');
m.addMaterial(materialdir,'copper');
%% Construction of Mesh
glib_sg_HalfSalientPole(5,8,2,10,4,0.5,(10:-1:0).^alpha/10^alpha,3*0.252,2*0.251,0.25);
m.read_g2d_bin('geom.g2d','MG0');
m.setMaterial('coil','copper');
m.setMaterial('rotor','m19_24ga');

m.ggmesh;

%% Solver
s = IHLTHSTL3(m);
s.units.length = 'mm';
s.units.currentDensity = 'A/mm^2';
s.units.magneticVectorPotential = 'A/m';
s.setExcitation('coil',100);
%%
% s.clearallbcs;
s.setdbc(s.m.getnIndexOnLine([0,0],[0,1]),1);
% calling solver

s.solve;

k = s.m.getnIndexOnLine([0,14.25],[1,14.25]);

p = s.m.nodes(k,:);
s.evalBn;
b = [s.Bn.x(k,:), s.Bn.y(k,:)];

[p(:,1),index] = sort(p(:,1));
p(:,2) = p(index,2);

b = b(index,:);
f = abs(fft(b(:,2)));
f = f(2:floor((length(f)-1)/2));

kr = sqrt(sum(f(2:end).^2))/f(1);

plot(p(:,1),b(:,2));
