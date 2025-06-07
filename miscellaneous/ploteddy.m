close all
t = s.m.t(:,1:3);
t = t';
t = t(:);
tt = 1:3*s.m.Ngt;
tt = reshape(tt,3,[]);
s = s.evalBxBy;
for nn  = 1:20
    for tk = 2:Nt+1
        Bk = sqrt(sum([s.Bx(:,tk),s.By(:,tk)].^2,2));
        Bk = repmat(Bk,1,3);
        Bk = Bk';
        Bk = Bk(:);
        
        
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
