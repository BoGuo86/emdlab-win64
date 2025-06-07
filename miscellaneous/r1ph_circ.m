

vs = @(t) sin(2*pi*50*t);
r = 1;

d1 = sc_diod;
d2 = sc_diod;
d3 = sc_diod;
d4 = sc_diod;


ts = 1e-5;
t = 0:ts:0.1;
Nt = length(t);

vo = zeros(1,Nt);

xu = 0;
xd = 0;
for i = 1:Nt-1
  

  r1 = d1.getR(vs(t(i)) - xu);
r2 = d2.getR(- xu);
r3 = d3.getR(xd - vs(t(i)));
r4 = d4.getR(xd);
  
  A = [1/r+1/r2+1/r1,-1/r;-1/r,1/r+1/r4+1/r3];
  
  x = A\[vs(t(i+1))/r1;vs(t(i+1))/r3];
  xu = x(1);
  xd = x(2);
  
  vo(i+1) = xu-xd;
end

plot(t,vo)


