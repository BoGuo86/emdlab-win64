clc
clear

vs = @(t) sin(2*pi*50*t);
r = 1;
l = 5e1;

ts = 1e-3;
t = 0:ts:0.08;
Nt = length(t);

p = tf_pulse;
p.t0 = 0;
p.duty = 0.1;

ty = sc_thyristor;

vo = zeros(1,Nt);
io = zeros(1,Nt);

for i = 1:Nt-1
  
  vt = vs(t(i)) - vo(i);
  rt = ty.getR(vt, p.getValue(t(i)), io(i));
  A = [1,rt;ts/2,-r*ts/2-l];
  x = A\[vs(t(i+1));-ts/2*vo(i)+r*(ts/2)*io(i)-l*io(i)];
  
  if x(2)<=0
%     vo(i) = 0;
%     io(i) = 0;
    x(1) = 0;
    x(2) = 0;
  end
  
  vo(i+1) = x(1);
  io(i+1) = x(2);
  
end


subplot(211)
plot(t,vo)
subplot(212)
plot(t,io)
