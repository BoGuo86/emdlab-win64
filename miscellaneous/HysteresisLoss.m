kh = 33.906897e-3;

alpha = 1;
beta = 2;

tetai = linspace(0,2*pi,1000);

k1 = kh/((2*pi)^(alpha-1))/(sum((abs(cos(tetai))).^alpha .* ...
    (abs(sin(tetai))).^(beta-alpha))*tetai(2));

bx = abs(Bxs);
by = abs(Bys);
dt = ((pi/4)/length(xrotorAngle))/w;
period = (pi/4)/w;

dBxdt = abs(Bxs(:,2:end)-Bxs(:,1:end-1))/dt;
dBydt = abs(Bys(:,2:end)-Bys(:,1:end-1))/dt;

phystx = (dBxdt.^alpha) .* (bx(:,1:end-1).^(beta-alpha));
phystx = sum(phystx ,2);
physty = (dBydt.^alpha) .* (by(:,1:end-1).^(beta-alpha));
physty = sum(physty ,2);
physt = phystx+physty;



TotalHyst = diag(s.m.gta(s.m.tzi(:,13))) * physt*7650;

TotalHyst = k1 * sum(TotalHyst) * 55e-9 * 4 * dt/period

physt = physt/4;

close all
t = s.m.t(s.m.tzi(:,13),1:3);
tt = 1:3*size(t,1);
t = t';
t = t(:);
tt = reshape(tt,3,[]);

close all
hold all
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

