%% Calculation of Hysteresis Loss Using MSE Method
% Bx, By and zoneIndex
angularSpeed = rpm * 2*pi/60;
% hysteresis loss coefficient
Kh = 33.906897e-3;
% power of frequency
alpha = 1;
% power of peak of magnetic flux density
beta = 2;
% calculation of KGSE
temp = linspace(0,2*pi,1000);
KGSE = Kh/((2*pi)^(alpha-1))/(sum((abs(cos(temp))).^alpha .* ...
    (abs(sin(temp))).^(beta-alpha))*temp(2));

bx = abs(Bxs(:,2:end-1));
by = abs(Bys(:,2:end-1));

dBxdt = abs(Bxs(:,3:end)-Bxs(:,1:end-2))/2;
dBydt = abs(Bys(:,3:end)-Bys(:,1:end-2))/2;

phystx = (dBxdt.^alpha) .* (bx.^(beta-alpha));
phystx = sum(phystx ,2);
physty = (dBydt.^alpha) .* (by.^(beta-alpha));
physty = sum(physty ,2);
physt = phystx+physty;

physt = physt*7650*KGSE* angularSpeed^alpha * length(xrotorAngle)^(alpha-1)/(simulationAngle)^alpha;

TotalHyst = diag(s.m.gta(s.m.tzi(:,13))) * physt;

TotalHyst = sum(TotalHyst) * 55e-9 * 4 

close all
hold all
t = s.m.t(s.m.tzi(:,13),1:3);
tt = 1:3*size(t,1);
t = t';
t = t(:);
tt = reshape(tt,3,[]);
c = colorbar;
c.Color = 'w';
c.Ticks = linspace(min(physt),max(physt),10);
c.Limits = [min(physt),max(physt)];
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

