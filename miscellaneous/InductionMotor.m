%% Simulation of induction motor
% this m file is for simulation of induction motor under free running
% condition, formulation is in stator reference frame

function InductionMotor

%% Initialization
clc;

%% Inputs
% source parameters
% line to line voltage [Volt]
Vll = 380;
% phase connection: star connection
Vph = Vll/sqrt(3);
% input frequency [Hz]
f = 60;
% base speed
wb = 2*pi*f;
% number of poles
p = 4;

% Machine parameters
% stator resistance
rs = 0.531;
% rotor resistance
rr = 0.408;
% mutual reactance
Xm = 31.95;
% lekage stator reactance
Xls = 0.95;
% lekage rotor reactance
Xlr = 0.95;
% stator reactance
Xss = Xls + Xm;
% rotor reactance
Xrr = Xlr + Xm;
D = Xss*Xrr - Xm^2;
% rotor inertia
J = 0.0089;
% friction coefficient
Damper = 0.001;

%% Simulation time
% step time: lower than 1e-4 give good results
ts = 1e-4;
% total simulation time
t = 0:ts:0.5;
% number of steps
Nt = length(t);

% input voltage
vs = Vph * [sin(2*pi*f*t); sin(2*pi*f*t - 2*pi/3); sin(2*pi*f*t - 4*pi/3)];

% qd0 stator voltage in stator reference frame
vsqd = getKs(0)*vs;

% allocation of memory for some variables
% qd0 stator currents
isqd = zeros(3, Nt);
% qd0 stator linkage flux
psisqd = zeros(3, Nt);
% qd0 rotor currents
irqd = zeros(3, Nt);
% qd0 rotor linkage flux
psirqd = zeros(3, Nt);
% qd0 rotor voltages
vrqd = zeros(3, Nt);

% electrical torque
Te = zeros(1, Nt);
% angular rotor speed
wr = zeros(1, Nt);

% a matrix that transform current to linkage flux
Psi2I = (1/D)*[Xrr, 0, 0, -Xm, 0, 0;
    0, Xrr, 0, 0, -Xm, 0;
    0, 0, D/Xls, 0, 0, 0;
    -Xm, 0, 0, Xss, 0, 0;
    0, -Xm, 0, 0, Xss, 0;
    0, 0, 0, 0, 0, D/Xlr];

%% Loop for calculation

for i = 2:Nt
    
    % first step: evaluation of linkage flux via forward integration
    psi_tmp = getPsi2V(0, wr(i-1))\([vsqd(:,i); vrqd(:,i)] + ([psisqd(:,i-1); psirqd(:,i-1)]/wb/ts));
    
    % second step: evaluation of currents
    i_tmp = Psi2I * psi_tmp;
    
    % putting data into variables
    psisqd(:,i) = psi_tmp(1:3);
    psirqd(:,i) = psi_tmp(4:6);
    
    isqd(:,i) = i_tmp(1:3);
    irqd(:,i) = i_tmp(4:6);
    
    % third step: evaluation of electrical torque
    Te(i) = (3/2)*(p/2)*(1/wb)*(psirqd(1,i)*irqd(2,i) - psirqd(2,i)*irqd(1,i));
    
    % forth step: solving mechanical equation
    wr(i) = (Te(i)*2/p + J*wr(i-1)/ts) / (J/ts + Damper); 
    
end

% inverse transform for evaluation of stator currents
is = getKsInverse(0)*isqd;

%% Plotting data
subplot(311)
plot(t, is)
title('Phase currents')

subplot(312)
plot(t, Te)
title('Electrical torque')

subplot(313)
plot(t, wr)
title('Rotor angular speed')

%% Nested functions

    function y = getPsi2V(omega_g, omega_r)
        y = [rs*Xrr/D + 1/wb/ts, omega_g/wb, 0, -rs*Xm/D, 0, 0;
            -omega_g/wb, rs*Xrr/D + 1/wb/ts, 0, 0, -rs*Xm/D, 0;
            0, 0, rs/Xls+1/wb/ts, 0, 0, 0;
            -rr*Xm/D, 0, 0, rr*Xss/D+1/wb/ts, (omega_g-omega_r)/wb, 0;
            0, -rr*Xm/D, 0, -(omega_g-omega_r)/wb,rr*Xss/D+1/wb/ts, 0;
            0, 0, 0, 0, 0, rr/Xlr+1/wb/ts];
    end

end

%% Auxiliary functions

function y = getKs(teta)
y = (2/3)*[cos(teta), cos(teta - 2*pi/3), cos(teta - 4*pi/3);
    sin(teta), sin(teta - 2*pi/3), sin(teta - 4*pi/3);
    1/2, 1/2, 1/2];
end

function y = getKsInverse(teta)
y = [cos(teta), sin(teta), 1;
    cos(teta - 2*pi/3), sin(teta - 2*pi/3), 1;
    cos(teta - 4*pi/3), sin(teta - 4*pi/3), 1];
end


