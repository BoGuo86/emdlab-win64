
function y = pulseFcn(t)

y = rem(t,0.02);
if y<0.01
  y=1;
else
  y = -1;
end

end
