function mz = MinimalMesh1(f,v)

check_g2d(f,v)
tmp = v(f(:,2),:)-v(f(:,1),:);
tmp = sqrt(sum(tmp.^2, 2));
maxLength = 1.05*max(tmp);

dt = delaunayTriangulation(v,f);
mz = TMZPC(dt.ConnectivityList(ext_dpoly2d(dt.incenter,f,v)<-1e-6,:),dt.Points);

aaa = 0;
while true
    ff = mz.edges(mz.bedges, :);
    tmp = mz.getCenterOfElements;
    tmp = tmp(mz.getMaxEdgeLength>maxLength | mz.getAspectRatio>3, :);
    if ~any(tmp)
        break
    end
    vv = [mz.nodes; tmp];
    dt = delaunayTriangulation(vv,ff);
    mz = TMZPC(dt.ConnectivityList(ext_dpoly2d(dt.incenter,f,v)<-1e-6,:),dt.Points);
    mz.moveNodes;
    aaa = aaa + 1;
    if aaa>5
        break
    end
end

end
