%% Calculation of Hysteresis Loss Using GSE Method
% Bx, By and zoneIndex
angularSpeed = p.rpm * 2*pi/60;
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
ZoneIndex = szi;
bx = Bxs;
by = Bys;

dbxdt = [bx(:,2)-bx(:,1),...
    (bx(:,3:end)-bx(:,1:end-2))/2,...
    bx(:,end)-bx(:,end-1)];

dbydt = [by(:,2)-by(:,1),...
    (by(:,3:end)-by(:,1:end-2))/2,...
    by(:,end)-by(:,end-1)];

phystx = (abs(dbxdt).^alpha) .* (abs(bx).^(beta-alpha));
physty = (abs(dbydt).^alpha) .* (abs(by).^(beta-alpha));

phystx(:,1) = phystx(:,1)/2;
phystx(:,end) = phystx(:,end)/2;

physty(:,1) = physty(:,1)/2;
physty(:,end) = physty(:,end)/2;

phystx = sum(phystx ,2);
physty = sum(physty ,2);
physt = phystx+physty;

physt = physt*7650*KGSE* angularSpeed^alpha * length(xrotorAngle)^(alpha-1)/(simulationAngle)^alpha;

TotalHyst = diag(s.m.gta(s.m.tzi(:,ZoneIndex))) * physt;

TotalHyst = sum(TotalHyst) * 55e-9 *4

close all
hold all
t = s.m.t(s.m.tzi(:,ZoneIndex),1:3);
tt = 1:3*size(t,1);
t = t';
t = t(:);
tt = reshape(tt,3,[]);
c = colorbar;
c.Color = 'w';
c.Ticks = linspace(min(physt),max(physt),8);
c.Limits = [min(physt),max(physt)];
c.FontSize = 20;
physt = repmat(physt,1,3);
physt = physt';
trisurf(tt',s.m.p(t,1),...
    s.m.p(t,2),physt(:),'edgecolor','none');
axis off equal;
view([0,0,1]);
colormap jet;
set(gcf,'Color','k')
title('Hysteresis Loss Density [W/m^3] (GSE Method)',...
    'color','w','fontsize',15);
set(gcf,'Renderer','opengl');
zoom on;

