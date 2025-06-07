%% File description:
% Short circuit force distribution
%% initialization
clc;
clear;
% path of folder containing materials data
materialdir = [cd,'\emdlab_release\Materials'];


w1 = 30;
w2 = 65; 
w3 = 30;
w4 = 65;
w5 = 30;

h1 = 109;
h2 = 1794;

meshSize = 15;

AT = 265581;

glib_tr1ph_ir0(w1,w2,w3,w4,w5,h1,h2,meshSize);

m = TMDBC;

% adding material to mesh data base
m.addMaterial(materialdir,'air');

m.read_g2d_bin('geom.g2d', 'mg3');


%% Solver
% getting an instance of solver object
% IHLTHSTL3: Isotropic Homogenous Linear THermo-Static Triangular Lagrangian 3 node
% each solver gets its compatible mesh as input
s = IHNLNRMSTL3(m);
% setting physical units
s.setUnit('length', 'mm');
s.setUnit('currentDensity', 'A/mm^2');
s.setUnit('magneticVectorPotential', 'A/m');


s.setExcitation('LV_coil', 'value', AT);
s.setExcitation('HV_coil', 'value', -AT);

m.ggmesh;
s.setAzBC(m.getfbniol_p0u([w1+w2+w3+w4+w5,0],[0,1]),0);

s.solve;

close all;

% s.plotBmagSmooth
s.evalBe
Fx = -s.edata.InternalCurrentDensity' .* s.results.Bey;
Fy = s.edata.InternalCurrentDensity' .* s.results.Bex;

m.showvf(Fx,Fy)
%%
clc

Fra = (4*pi*1e-7/2)*AT^2/h2/1e-3;
fprintf('Fr (analytical) = %f\n', Fra);

mz = m.mzs.LV_coil;

tmp = sum(Fx(m.ezi(:,mz.zi)).*mz.getAreaOfElements);
fprintf('Fr_LV (numerical) = %f\n', tmp);
tmp = sum(abs(Fy(m.ezi(:,mz.zi))).*mz.getAreaOfElements)/2;


mz = m.mzs.HV_coil;

tmp = sum(Fx(m.ezi(:,mz.zi)).*mz.getAreaOfElements);
fprintf('Fr_HV (numerical) = %f\n', tmp);
tmp = sum(abs(Fy(m.ezi(:,mz.zi))).*mz.getAreaOfElements)/2;




% F = s.edata.InternalCurrentDensity .* repmat((m.gea / 3), 3, 1);
% % applying scales on load vector
% F = sparse(m.cl', ones(3 * m.Ne, 1), F);
% F= full(F);
% s.evalBn;
% Fx = -F .* s.results.Bny;
% Fy = F .* s.results.Bnx;


% fprintf('Fa_LV (numerical) = %f\n', sum(abs(Fy(m.mzs.LV_coil.l2g)))/2);
% fprintf('Fa_HV (numerical) = %f\n', sum(abs(Fx(m.mzs.HV_coil.l2g)))/2);