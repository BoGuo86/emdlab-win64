clc

Rph = 3;
D = 0.01;
Fload = 100;
Lm = 12e-3;
mass = 2;

ts = 1e-4;
tf = 0.4;
t = 0:ts:tf;
Nt = length(t);

A = [-Rph/Lm,0,-kem/Lm;0,0,1;2*kem/mass,0,-D/mass];
B = [1/2/Lm,0;0,0;0,-1/mass];

myfun = @(t,y) A*y+B*[getVdc_step(t);Fload];
[t,y] = ode45(myfun, t, [0;0;0]);

plot(t, v)
box on;
title('Speed [m/s]')
hold on
plot(t,y(:,3))


function y = getVdc_step(t)

    if t > 0.1 && t < 0.15
        y = 55;
    elseif t > 0.15 && t < 0.2
        y = 40;
    else
        y = 48;
    end

end
