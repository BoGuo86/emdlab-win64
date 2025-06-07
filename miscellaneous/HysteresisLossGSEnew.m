
function [sloss,rloss] = HysteresisLossGSEnew(Bx,By,s,xi)

%% Calculation of Hysteresis Loss Using GSE Method
% Bx, By and zoneIndex
rpm = 1500;
simulationAngle = pi/2;
angularSpeed = rpm * 2*pi/60;
% hysteresis loss coefficient
Kh =  6.673865e-3;
% power of frequency
alpha = 1.2916;
% power of peak of magnetic flux density
beta = 2;
% calculation of KGSE
temp = linspace(0,2*pi,1000);
KGSE = Kh/((2*pi)^(alpha-1))/(sum((abs(cos(temp))).^alpha .* ...
    (abs(sin(temp))).^(beta-alpha))*temp(2));

% zone
szi = 10;
ZoneIndex = szi;
bx = Bx(s.m.tzi(:,ZoneIndex),:);
by = By(s.m.tzi(:,ZoneIndex),:);

ntmp = length(xi);
myint = diff(xi);

tmmp = myint(1:end-1)+myint(2:end);
tmmp = sparse(1:ntmp-2,1:ntmp-2,1./tmmp);
dbxdt = [diff(bx(:,1:2),1,2)/myint(1),(-bx(:,1:end-2)+bx(:,3:end))*tmmp,...
    diff(bx(:,end-1:end),1,2)/myint(end)];

dbydt = [diff(by(:,1:2),1,2)/myint(1),(-by(:,1:end-2)+by(:,3:end))*tmmp,...
    diff(by(:,end-1:end),1,2)/myint(end)];

phystx = (abs(dbxdt).^alpha) .* (abs(bx).^(beta-alpha));
physty = (abs(dbydt).^alpha) .* (abs(by).^(beta-alpha));

phystx = (phystx(:,1:end-1)+phystx(:,2:end))/2;
physty = (physty(:,1:end-1)+physty(:,2:end))/2;

phystx = phystx*sparse(1:ntmp-1,1:ntmp-1,myint);
physty = physty*sparse(1:ntmp-1,1:ntmp-1,myint);

phystx = sum(phystx ,2);
physty = sum(physty ,2);
physt = phystx+physty;


% physt = physt*7650*KGSE* angularSpeed^alpha * length(xi)^(alpha-1)/(simulationAngle)^alpha;
physt = physt*7650*KGSE* angularSpeed^alpha/simulationAngle;

sloss = diag(s.m.gta(s.m.tzi(:,ZoneIndex))) * physt;

sloss = sum(sloss) * 55e-9 *4;

szi = 9;
ZoneIndex = szi;
bx = Bx(s.m.tzi(:,ZoneIndex),:);
by = By(s.m.tzi(:,ZoneIndex),:);

ntmp = length(xi);
myint = diff(xi);

tmmp = myint(1:end-1)+myint(2:end);
tmmp = sparse(1:ntmp-2,1:ntmp-2,1./tmmp);
dbxdt = [diff(bx(:,1:2),1,2)/myint(1),(-bx(:,1:end-2)+bx(:,3:end))*tmmp,...
    diff(bx(:,end-1:end),1,2)/myint(end)];

dbydt = [diff(by(:,1:2),1,2)/myint(1),(-by(:,1:end-2)+by(:,3:end))*tmmp,...
    diff(by(:,end-1:end),1,2)/myint(end)];

phystx = (abs(dbxdt).^alpha) .* (abs(bx).^(beta-alpha));
physty = (abs(dbydt).^alpha) .* (abs(by).^(beta-alpha));

phystx = (phystx(:,1:end-1)+phystx(:,2:end))/2;
physty = (physty(:,1:end-1)+physty(:,2:end))/2;

phystx = phystx*sparse(1:ntmp-1,1:ntmp-1,myint);
physty = physty*sparse(1:ntmp-1,1:ntmp-1,myint);

phystx = sum(phystx ,2);
physty = sum(physty ,2);
physt = phystx+physty;


% physt = physt*7650*KGSE* angularSpeed^alpha * length(xi)^(alpha-1)/(simulationAngle)^alpha;
physt = physt*7650*KGSE* angularSpeed^alpha/simulationAngle;

rloss = diag(s.m.gta(s.m.tzi(:,ZoneIndex))) * physt;

rloss = sum(rloss) * 55e-9 *4;
end

% close all
% hold all
% t = s.m.t(s.m.tzi(:,ZoneIndex),1:3);
% tt = 1:3*size(t,1);
% t = t';
% t = t(:);
% tt = reshape(tt,3,[]);
% c = colorbar;
% c.Color = 'w';
% c.Ticks = linspace(min(physt),max(physt),8);
% c.Limits = [min(physt),max(physt)];
% c.FontSize = 20;
% physt = repmat(physt,1,3);
% physt = physt';
% trisurf(tt',s.m.p(t,1),...
%     s.m.p(t,2),physt(:),'edgecolor','none');
% axis off equal;
% view([0,0,1]);
% colormap jet;
% set(gcf,'Color','k')
% title('Hysteresis Loss Density [W/m^3] (GSE Method)',...
%     'color','w','fontsize',15);
% set(gcf,'Renderer','opengl');
% zoom on;

