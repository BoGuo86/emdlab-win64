
timeStep = 0.01;
Ntime = 300;
L = 30;
Nx = 100;
x = linspace(0, L, Nx);
xStep = x(2);

U = zeros(Ntime, Nx);

k = 1;
q = 0.24;
v = 1.506;
xp1 = 20;
xp2 = 25;

for k = 1:Nx
  U(1,k) = (3/8)*(2-1/v)*sech(x(k)-xp1)^2 + sech(x(k)-xp2)^2;
end

U(1,1) = 0;
U(1,end) = 0;

for n = 1:Ntime-1
  
  for k = 2:Nx-1
    A = -2*q/timeStep^2;
    B = 1/timeStep + 4*q*U(n,k)/timeStep^2;
    C = -U(n,k)/timeStep - ...
      2*q*U(n,k)^2/timeStep^2 - ...
      k+x(*(U(n,k+1)-U(n,k))/xStep;
    
    delta = B^2 - 4*A*C;
    
    Utmp = (-B + sqrt(delta))/(2*A);
    if isreal(Utmp) && Utmp>0
      U(n+1,k) = Utmp;
    else
      Utmp = (-B - sqrt(delta))/(2*A);
      if isreal(Utmp) && Utmp>0
        U(n+1,k) = Utmp;
      end
    end
  end

end
[xgrid, tgrid] = meshgrid(x,0:timeStep:(Ntime-1)*timeStep);
surf(xgrid, tgrid, U, 'EdgeColor', 'None')
view([0,0,1])
