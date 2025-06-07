%% Calculation of Eddy Current Loss in Electrical Lamination
% needed parameters
angularSpeed = rpm * 2*pi/60;
% zone index
zoneIndex = 12;
% depth of motor
depth = 55e-3;

dBxdt = (Bxr(:,3:end)-Bxr(:,1:end-2))/2;
dBydt = (Byr(:,3:end)-Byr(:,1:end-2))/2;

peddy = (dBxdt.^2+dBydt.^2);
peddy = sum(peddy,2);
peddy =  (78.856601e-6/2/pi^2)*peddy*7650*length(xrotorAngle)*angularSpeed^2/(simulationAngle)^2;

TotalEddy = diag(s.m.gta(s.m.tzi(:,zoneIndex))) * peddy;

TotalEddy = sum(TotalEddy) * depth * 1e-6 * 4

close all
hold all

t = s.m.t(s.m.tzi(:,zoneIndex),1:3);
tt = 1:3*size(t,1);
t = t';
t = t(:);
tt = reshape(tt,3,[]);

c = colorbar;
c.Color = 'w';
c.Ticks = linspace(min(peddy),max(peddy),10);
c.Limits = [min(peddy),max(peddy)];

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