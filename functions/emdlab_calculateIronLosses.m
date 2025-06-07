function Pfe = emdlab_calculateIronLosses(mzptr, Bex, Bey, simTime, Lstk)

% density of lamination
Density = 7850;

% eddy current coefficient
Ke = 3.4505e-05;

% Hysteresis loss coefficient
Kh =  5.3941e-03 ;

% power of frequency
alpha = 1.3092e+00;

% power of magnetic field
beta = 2;

% calculation of KGSE
tmp = linspace(0,2*pi,1000);
KGSE = Kh/((2*pi)^(alpha-1))/(sum((abs(cos(tmp))).^alpha.*(abs(sin(tmp))).^(beta-alpha))*tmp(2));

% eddy current loss
Peddy = int_trap((diff_cen_mean(Bex,simTime).^2 + diff_cen_mean(Bey,simTime).^2), simTime);
Peddy = (Ke/2/pi^2) * (1/simTime(end)) * (Peddy.*mzptr.getAreaOfElements) * Lstk * 1e-9 * Density;

% hysteresis loss
Physt = int_trap(abs(diff_cen_mean(Bex,simTime)).^alpha.*abs(Bex).^(beta-alpha) ...
    + abs(diff_cen_mean(Bey,simTime)).^alpha.*abs(Bey).^(beta-alpha), simTime);
Physt = KGSE*(1/simTime(end))*(Physt.*mzptr.getAreaOfElements) * Lstk * 1e-9 * Density;

% total iron loss
Pfe = Peddy + Physt;

% clc
% disp('Total rotor eddy loss:');
% disp(sum(Peddy));
mzptr.smoothPlot(Pfe*1e9/Lstk./mzptr.getAreaOfElements/Density)

Pfe = sum(Peddy + Physt);

end