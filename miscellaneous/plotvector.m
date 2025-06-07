close all
tr = triangulation(s.m.t(:,1:3),s.m.p);
c = tr.incenter;
s = s.evalBxBy;
for nn  = 1:20
    for tk = 2:Nt+1
        quiver(c(:,1),c(:,2),s.Bx(:,tk),s.By(:,tk),1,'color','w')
        axis off equal;
        axis([20 50 20 50])
        view([0,0,1]);
        colormap hot;
        set(gcf,'Color','k')

        set(gcf,'Renderer','opengl');

        pause(0.05)
    end
end
