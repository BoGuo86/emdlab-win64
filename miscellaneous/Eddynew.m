%% Calculation of Eddy Current Loss in Electrical Lamination
% needed parameters
angularSpeed = 1500 * 2*pi/60;
% zone index
zoneIndex =szi;
bx = spline(xrotorAngle,Bxs);
by = spline(xrotorAngle,Bys);
% depth of motor
depth = 55e-3;

tmp = length(xrotorAngle);
cc = ones(1,tmp-1);
cc = [cc,0]+[0,cc];


dBxdt = bx;
dBxdt.coefs = dBxdt.coefs*diag(3:-1:1,1);
dBydt = by;
dBydt.coefs = dBydt.coefs*diag(3:-1:1,1);

peddy = (ppval(dBxdt,xrotorAngle).^2+ppval(dBydt,xrotorAngle).^2);
peddy = peddy*diag(cc);
peddyt = sum(peddy);
peddy = sum(peddy,2)*xrotorAngle(2)/2;
peddy =  (66.156534e-6/2/pi^2)*peddy*7650*angularSpeed^2/simulationAngle;

TotalEddy = diag(s.m.gta(s.m.tzi(:,zoneIndex))) * peddy;

TotalEddy = sum(TotalEddy) * depth * 1e-6 

close all
hold all

t = s.m.t(s.m.tzi(:,zoneIndex),1:3);
tt = 1:3*size(t,1);
t = t';
t = t(:);
tt = reshape(tt,3,[]);

c = colorbar;
c.Color = 'w';
c.Ticks = linspace(min(peddy),max(peddy),8);
c.Limits = [min(peddy),max(peddy)];
c.FontSize = 20;

peddy = repmat(peddy,1,3);
peddy = peddy';

trisurf(tt',s.m.p(t,1),...
    s.m.p(t,2),peddy(:),'edgecolor','none');

axis off equal;
view([0,0,1]);
colormap jet;
set(gcf,'Color','k')
title('Eddy Current Loss Density [W/m^3]',...
    'color','w','fontsize',15);
set(gcf,'Renderer','opengl');
zoom on;