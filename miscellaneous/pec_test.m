
clc
clear

tic
c = pec_db(3,3);

c.stickElement(1,1,2);
c.stickElement(2,2,3);
c.stickElement(3,3,1);

c.setElementObject(1, pelib_ivs_sinusoidal);
c.setElementObject(2, pelib_idiod);
c.setElementObject(3, pelib_resistance);

c.stickGround2Node(1);

ts = 1e-4;
tf = 0.06;
t = 0:ts:tf;
Nt = length(t);

vi = zeros(1,Nt);
vo = zeros(1,Nt);

for i = 2:Nt
    
end

toc