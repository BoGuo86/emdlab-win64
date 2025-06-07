%% Calculation of Eddy Current Loss in Electrical Lamination
% needed parameters
function [sloss,rloss] = EdyyCurrentNew(Bx,By,s,xi)
rpm = 1500;
angularSpeed = 1500 * 2*pi/60;
simulationAngle = pi/2;
% depth of motor
depth = 55e-3;
% zone index
zoneIndex = 10;
bx = Bx(s.m.tzi(:,zoneIndex),:);
by = By(s.m.tzi(:,zoneIndex),:);

ntmp = length(xi)-1;
myint = diff(xi);

dbxdt = diff(bx,1,2)*sparse(1:ntmp,1:ntmp,1./myint);
dbxdt = [dbxdt(:,1),(dbxdt(:,1:end-1)+dbxdt(:,2:end))/2,dbxdt(:,end)];

dbydt = diff(by,1,2)*sparse(1:ntmp,1:ntmp,1./myint);
dbydt = [dbydt(:,1),(dbydt(:,1:end-1)+dbydt(:,2:end))/2,dbydt(:,end)];

peddy = (dbxdt.^2+dbydt.^2);
peddy = (peddy(:,1:end-1)+peddy(:,2:end))/2;

peddy = peddy*sparse(1:ntmp,1:ntmp,myint);

% peddy = peddy(:,1)+2*sum(peddy(:,3:2:2*tmp-1))+4*sum(peddy(:,2:2:2*tmp))

peddy = sum(peddy,2);
peddy =  (66.156534e-6/2/pi^2)*peddy*7650*angularSpeed^2/simulationAngle;

TotalEddy = diag(s.m.gta(s.m.tzi(:,zoneIndex))) * peddy;

sloss = sum(TotalEddy) * depth * 1e-6 * 4;

zoneIndex = 9;
bx = Bx(s.m.tzi(:,zoneIndex),:);
by = By(s.m.tzi(:,zoneIndex),:);

ntmp = length(xi)-1;
myint = diff(xi);

dbxdt = diff(bx,1,2)*sparse(1:ntmp,1:ntmp,1./myint);
dbxdt = [dbxdt(:,1),(dbxdt(:,1:end-1)+dbxdt(:,2:end))/2,dbxdt(:,end)];

dbydt = diff(by,1,2)*sparse(1:ntmp,1:ntmp,1./myint);
dbydt = [dbydt(:,1),(dbydt(:,1:end-1)+dbydt(:,2:end))/2,dbydt(:,end)];

peddy = (dbxdt.^2+dbydt.^2);
peddy = (peddy(:,1:end-1)+peddy(:,2:end))/2;

peddy = peddy*sparse(1:ntmp,1:ntmp,myint);

% peddy = peddy(:,1)+2*sum(peddy(:,3:2:2*tmp-1))+4*sum(peddy(:,2:2:2*tmp))

peddy = sum(peddy,2);
peddy =  (66.156534e-6/2/pi^2)*peddy*7650*angularSpeed^2/simulationAngle;

TotalEddy = diag(s.m.gta(s.m.tzi(:,zoneIndex))) * peddy;

rloss = sum(TotalEddy) * depth * 1e-6 * 4;

% close all
% hold all

% t = s.m.t(s.m.tzi(:,zoneIndex),1:3);
% tt = 1:3*size(t,1);
% t = t';
% t = t(:);
% tt = reshape(tt,3,[]);
% 
% c = colorbar;
% c.Color = 'w';
% c.Ticks = linspace(min(peddy),max(peddy),8);
% c.Limits = [min(peddy),max(peddy)];
% c.FontSize = 20;
% 
% peddy = repmat(peddy,1,3);
% peddy = peddy';

% trisurf(tt',s.m.p(t,1),...
%     s.m.p(t,2),peddy(:),'edgecolor','none');
% 
% axis off equal;
% view([0,0,1]);
% colormap jet;
% set(gcf,'Color','k')
% title('Eddy Current Loss Density [W/m^3]',...
%     'color','w','fontsize',15);
% set(gcf,'Renderer','opengl');
% zoom on;

end
