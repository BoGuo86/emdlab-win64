%% setting x dimension
L = 20;
Nx = 200;
x = linspace(0, L, Nx);
xStep = x(2);
%% setting time
timeStep = 0.003;
Ntime = 101;
%% solution variable
U = zeros(Ntime, Nx);
%% inputs
kk = 10;
gamma = 0.001;
v = 1.506;
xp1 = 4;
xp2 = 12;
%% setting initial condition
for n = 1:2
  for k = 1:Nx
    U(n,k) = (3/8)*(2-1/v)*sech(x(k)-xp1)^2 + sech(x(k)-xp2)^2;
  end
end
%% setting boundary equation
U(1,1) = 0; U(2,1) = 0;
U(1,end) = 0; U(2,end) = 0;
%% loop: implicit method
A = zeros(Nx,Nx);
jIndex = 1;
for k = 2:Nx-1
  A(k, jIndex) = -gamma/xStep^2;
  A(k, jIndex+1) = 1+2*gamma/xStep^2;
  A(k, jIndex+2) = -gamma/xStep^2;
  jIndex = jIndex + 1;
end
% for derivation of boundary conditions we impose
A(1,1) = 1;
A(end,end) = 1;
b = zeros(Nx,1);
for n = 2:Ntime-1
  for k = 2:Nx-1 
    % Uxx: (n,k-1) <-> (n,k) <-> (n,k+1)
    Uxx = (U(n,k+1) - 2*U(n,k) + U(n,k-1))/xStep^2;
    % Ux: central
    Ux = (U(n,k+1) - U(n,k-1))/2/xStep;
    % Utb: first order backward
    Utb = (U(n,k) - U(n-1,k))/timeStep;
    % derivation in new time
    b(k) = U(n,k) + gamma*Uxx + timeStep*(-U(n,k) + Uxx + kk*Ux);
  end
  U(n+1,:) = A\b;
end
[xgrid, tgrid] = meshgrid(x,0:timeStep:(Ntime-1)*timeStep);
surf(xgrid, tgrid, U, 'EdgeColor', 'None')
colormap jet
view([0,0,1])
xlabel('x')
ylabel('t')
title('U + U_t = U_{xx} + \gamma*U_{txx} + k*U_x')
