function mzptr = emdlab_m2d_getPolygonTM(p)

f = (1:size(p,1))';
f = [f,circshift(f,-1)];
mzptr = emdlab_m2d_mm(f, p);

end