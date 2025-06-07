alpha = 1;
%% initialization
materialdir = [cd,'\MaterialsData'];
%% Mesh
m = TMDBC;
m.addMaterial(materialdir,'air');
m.addMaterial(materialdir,'m19_24ga');
m.addMaterial(materialdir,'copper');
%% Construction of Mesh
%         tmp = (10:-0.5:0).^alpha/10^alpha;
xx = 0:0.1:8;
%         % % tmp = exp(-alpha*xx)-exp(-alpha*8);
tmp = alpha*(1./cos(xx*pi/24)-1);
tmp = fliplr(tmp);

glib_sg_HalfSalientPole(5,10,2,10,4,0.25,tmp,1,0.6,0.25);
m.read_g2d_bin('geom.g2d','MG1',20);
m.setMaterial('coil','copper');
m.setMaterial('rotor','m19_24ga');
m.ggmesh;
%% Solver
s = IHNLNRMSTL3(m);
s.setUnit('length', 'mm');
s.setUnit('currentDensity', 'A/mm^2');
s.setUnit('magneticVectorPotential', 'A/m');
s.setExcitation('coil','value',-100,'type','c');
%% calling solver
s.bcs.clearAllBCs;
s.bcs.setDirichlet(s.m.getnIndexOnLine([0,0],[0,1]),0);
s.solve;
%% Fourier
k = s.m.getnIndexOnLine([0,14.25],[1,14.25]);
p = s.m.nodes(k,:);
s.evalBn;
b = s.results.Bny(k,:);
[~,index] = sort(p(:,1));
b = b(index)';
b = [b(1:end-1),-fliplr(b)];
b = [b(1:end-1),fliplr(b)];
h = abs(fft(b));
h = h(2:floor((length(h)-1)/2));