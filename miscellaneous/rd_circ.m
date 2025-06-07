

function rd_circ()

vs = @(t) sin(2*pi*50*t);
r = 1;

ron = 0.01;
roff = 1e6;

t = linspace(0,0.04,500);
Nt = length(t);

vo = zeros(1,Nt);

for i = 1:Nt-1
  
  vd = vs(t(i)) - vo(i);
  if vd>0
    vo(i+1) = vs(t(i+1)) * r/(r+ron);
  elseif vd<0
    vo(i+1) = vs(t(i+1)) * r/(r+roff);
  else
    vo(i+1) = vo(i);
  end
  
end

plot(t,vo)

