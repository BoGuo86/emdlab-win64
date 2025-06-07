close all
index = s.m.tzi(:,[18 19]);
index = logical(sum(index,2));
Bk = sqrt(sum((s.B(index,:)).^2,2));
Bk = repmat(Bk,1,3);
Bk = Bk';
Bk = Bk(:);
t = s.m.t(index,1:3);
t = t';
t = t(:);
tt = 1:3*sum(index);
tt = reshape(tt,3,[]);
trisurf(tt',s.m.p(t,1),...
    s.m.p(t,2),Bk,'edgecolor','none');
axis off equal;
view([0,0,1]);
c = colorbar;
c.Color = 'w';
c.Ticks = linspace(min(Bk),max(Bk),10);
c.Limits = [min(Bk),max(Bk)];
colormap hot;
set(gcf,'Color','k')
title('Magnetic Flux Density Amplitude [Tesla]',...
    'color','w','fontsize',15);
set(gcf,'Renderer','opengl');
shading interp
zoom on;