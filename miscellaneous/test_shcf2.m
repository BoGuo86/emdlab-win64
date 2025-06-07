%% COMMAN DATA
clc
clear

mu0 = 4*pi*1e-7;
DS = 1090;
HF = 2248;
ES = 2019;
ESR = 0;
BS = 1000;
SS = 923;
BJ = 1000;

uk_ind = 14.64;
uk_ohm = 0.17;

Sigma = 100;
E_Modul = 113000;

Ustand = 0.9;

%% CALCULATION OF FORCE OF EACH PARTS 

AT = 265581;
Di = 1150;
Bw = 65;
Aum = 109;
L0 = 90;
Hwm = 1794;
Bot = 2.2;
Hot = 8.3;

Dm = Di+2*Bw/3;
Frp = mu0*AT^2*pi*Dm/Hwm/2;

po = Frp/pi/Dm;

ri = Di/2;
ro = ri+Bw;

tmp = -po*ro^2/(ro^2-ri^2)*(ro^2/ri^2+1);

disp(tmp)
