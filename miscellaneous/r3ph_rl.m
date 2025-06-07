

va = @(t) 310*sin(2*pi*50*t);
vb = @(t) 310*sin(2*pi*50*t-2*pi/3);
vc = @(t) 310*sin(2*pi*50*t-4*pi/3);

r = 10;
l = 5e-2;

d1 = sc_Diod;
d2 = sc_Diod;
d3 = sc_Diod;
d4 = sc_Diod;
d5 = sc_Diod;
d6 = sc_Diod;

ts = 1e-5;
t = 0:ts:0.08;
Nt = length(t);

vo = zeros(1,Nt);
io = zeros(1,Nt);
xu = 0;
xd = 0;

for i = 1:Nt-1
  
  r1 = d1.getR(va(t(i)) - xu);
  r2 = d2.getR(vb(t(i)) - xu);
  r3 = d3.getR(vc(t(i)) - xu);
  r4 = d4.getR(xd - va(t(i)));
  r5 = d5.getR(xd - vb(t(i)));
  r6 = d6.getR(xd - vc(t(i)));
  
  A = [1/r1+1/r2+1/r3,0,1;0,1/r5+1/r4+1/r6,-1;ts/2,-ts/2,-r*ts/2-l];
  
  x = A\[va(t(i+1))/r1+vb(t(i+1))/r2+vc(t(i+1))/r3;...
    va(t(i+1))/r4+vb(t(i+1))/r5+vc(t(i+1))/r6;...
    -ts/2*(xu-xd)+r*ts/2*io(i)-l*io(i)];
  
  xu = x(1);
  xd = x(2);
  
  vo(i+1) = xu-xd;
  io(i+1) = x(3);
  
end

plot(t,io)


