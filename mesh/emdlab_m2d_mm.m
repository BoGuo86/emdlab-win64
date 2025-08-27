function mz = emdlab_m2d_mm(f, v, ces)

emdlab_g2d_validateFV(f,v)

if nargin == 2
    dt = delaunayTriangulation(v, f);
elseif nargin == 3
    dt = delaunayTriangulation(v, [f; ces]);
end

mz = emdlab_m2d_tmz(dt.ConnectivityList(ext_dpoly2d(dt.incenter,f,v)<-1e-6,:), dt.Points);

end
