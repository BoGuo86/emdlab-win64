%% setting x dimension
L = 30;
Nx = 120;
x = linspace(0, L, Nx);
xStep = x(2);
%% setting time
timeStep = 0.01;
Ntime = 301;
%% solution variable
U = zeros(Ntime, Nx);
%% inputs
kk = 1;
q = 0.24;
v = 1.506;
xp1 = 20;
xp2 = 25;
%% setting initial condition
for n = 1:2
  for k = 1:Nx
    U(n,k) = (3/8)*(2-1/v)*sech(x(k)-xp1)^2 + sech(x(k)-xp2)^2;
  end
end
%% setting boundary equation
U(1,1) = 0; U(2,1) = 0;
U(1,end) = 0; U(2,end) = 0;
%% loop: explicit method
for n = 2:Ntime-1
  for k = 2:Nx-1
    % Uxx: (n,k-1) <-> (n,k) <-> (n,k+1)
    Uxx = (U(n,k+1) - 2*U(n,k) + U(n,k-1))/xStep^2;
    % Ux: central
    Ux = (U(n,k+1) - U(n,k-1))/2/xStep;
    % Utb: first order backward
    Utb = (U(n,k) - U(n-1,k))/timeStep;
    % derivation in new time
    U(n+1,k) = U(n,k) + timeStep*(-U(n,k) + Uxx + 2*q*Utb^2 + kk*Ux);
  end
end
[xgrid, tgrid] = meshgrid(x,0:timeStep:(Ntime-1)*timeStep);
surf(xgrid, tgrid, U, 'EdgeColor', 'None')
colormap jet
view([0,0,1])
xlabel('x')
ylabel('t')
title('U + U_t = U_{xx} + 2*q*U_t^2 + k*U_x')
