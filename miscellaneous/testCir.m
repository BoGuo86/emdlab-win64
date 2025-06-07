

ts = 1e-4;
t = 0:ts:0.08;
Nt = length(t);

p = tf_pulse;
p.t0 = 0.005;
p.duty = 0.1;

y = zeros(1,Nt);

for i = 1:Nt
  y(i) = p.getValue(t(i));
end

plot(t,y)
