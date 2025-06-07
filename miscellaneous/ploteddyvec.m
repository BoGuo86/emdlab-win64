close all
scale = 1;

for nn  = 1:20
    for tk = 1:Nt
        quiver3(s.m.p(:,1),s.m.p(:,2),zeros(s.m.Ngp,1),...
                zeros(s.m.Ngp,1),zeros(s.m.Ngp,1),...
                s.ee(:,tk),scale,'color','w',...
                'MaxHeadSize',0.05);
        
        axis off equal;
        view([1,1,1]);
        colormap hot;
        set(gcf,'Color','k')

        set(gcf,'Renderer','opengl');

        pause(0.05)
    end
end