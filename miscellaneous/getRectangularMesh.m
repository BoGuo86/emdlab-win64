function mz = getRectangularMesh(x1,y1,x2,y2,meshSize)

Nx = ceil(abs(x2-x1)/meshSize);
Ny = ceil(abs(y2-y1)/meshSize);

x = linspace(x1,x2,Nx);
y = linspace(y1,y2,Nx);

[x,y] = meshgrid(x,y);

p = [x(:),y(:)];
t = delaunay(p);

mz = mlib_tmz(t,p);

end