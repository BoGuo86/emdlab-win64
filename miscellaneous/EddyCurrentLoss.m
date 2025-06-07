%% Calculation of Eddy Current Loss in Electrical Lamination
% needed parameters
angularSpeed = p.rpm * 2*pi/60;
% depth of motor
depth = 55e-3;
% zone index
zoneIndex = szi;
bx = Bxs;
by = Bys;

dbxdt = [bx(:,2)-bx(:,1),...
    (bx(:,3:end)-bx(:,1:end-2))/2,...
    bx(:,end)-bx(:,end-1)];

dbydt = [by(:,2)-by(:,1),...
    (by(:,3:end)-by(:,1:end-2))/2,...
    by(:,end)-by(:,end-1)];

peddy = (dbxdt.^2+dbydt.^2);

tmp = length(xrotorAngle);
% peddy = peddy(:,1)+2*sum(peddy(:,3:2:2*tmp-1))+4*sum(peddy(:,2:2:2*tmp))
peddy(:,1) = peddy(:,1) / 2;
peddy(:,end) = peddy(:,end) / 2;

peddy = sum(peddy,2);
peddy =  (66.156534e-6/2/pi^2)*peddy*7650*length(xrotorAngle)*angularSpeed^2/(simulationAngle)^2;

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