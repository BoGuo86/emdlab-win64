

ts = 1e-5;
t = 0:ts:0.08;
Nt = length(t);

vu = 5;
vd = -5;

r=1;

s1 = sc_md;
s2 = sc_md;
s3 = sc_md;
s4 = sc_md;
s5 = sc_md;
s6 = sc_md;

v = zeros(4,Nt);


for i = 1:Nt-1
  

  c = areFcn(t(i)) - sin(2*pi*50*t(i));
  vtmp = vu - v(1,i);
  r1 = s1.getR(c,vtmp);
  vtmp = v(1,i) - vd;
  r2 = s2.getR(-c,vtmp);
  
  c = areFcn(t(i)) - sin(2*pi*50*t(i)-2*pi/3);
  vtmp = vu - v(2,i);
  r3 = s3.getR(c,vtmp);
  vtmp = v(2,i) - vd;
  r4 = s4.getR(-c,vtmp);
  
  c = areFcn(t(i)) - sin(2*pi*50*t(i)-4*pi/3);
  vtmp = vu - v(3,i);
  r5 = s5.getR(c,vtmp);
  vtmp = v(3,i) - vd;
  r6 = s6.getR(-c,vtmp);
  
  A = [1/r+1/r1+1/r2,0,0,-1/r;...
    0,1/r+1/r3+1/r4,0,-1/r;...
    0,0,1/r+1/r5+1/r6,-1/r;...
    1,1,1,-3];
  
  b = [vu/r1+vd/r2;vu/r3+vd/r4;vu/r5+vd/r6;0];
  
  v(:,i+1) = A\b;
  
end

subplot(411)
plot(t,v(1,:))
subplot(412)
plot(t,v(2,:))
subplot(413)
plot(t,v(3,:))
subplot(414)
plot(t,v(4,:))
