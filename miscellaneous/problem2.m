
function problem2

x = -10:0.01:10;
plot(x,x.*cos(x),x,ones(1,length(x)));
title('Selected interval = [4,6]');

f = @(x) x*cos(x)-1;
bisection(f,4,6)


f = @(x) 1/cos(x);
fixedPoint(f,5)

function p = bisection(f,a,b)

if f(a)*f(b)>0
  disp('Wrong choice bro')
else
  p = (a + b)/2;
  err = abs(f(p));
  while err > 1e-7
    if f(a)*f(p)<0
      b = p;
    else
      a = p;
    end
    p = (a + b)/2;
    err = abs(f(p));
  end
end

function xR = fixedPoint(f,xR0)

xR = f(xR0);
err = abs(xR-xR0);

while err>1e-7
  xR0 = xR;
  xR = f(xR0);
  err = abs(xR-xR0);
end
