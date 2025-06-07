close all
t = s.m.t(:,1:3);
t = t';
t = t(:);
tt = 1:3*s.m.Ngt;
tt = reshape(tt,3,[]);
for nn  = 1:20
    for tk = 1:Nt
        Bk = sqrt(sum((s.EC(:,2*tk-1:2*tk)).^2,2));
        Bk = repmat(Bk,1,3);
        Bk = Bk';
        Bk = Bk(:);
        Bk = Bk/max(Bk);
        
        trisurf(tt',s.m.p(t,1),...
            s.m.p(t,2),Bk,'edgecolor','none');
        axis off equal;
        view([0,0,1]);
        colormap hot;
        set(gcf,'Color','k')

        set(gcf,'Renderer','opengl');

        pause(0.05)
    end
end