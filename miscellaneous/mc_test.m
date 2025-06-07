
clc
clear

tmp = linspace(0,2*pi,100);
tmp = tmp(1:end-1)';


m = mc_circularAirGap([20*cos(tmp)+1,20*sin(tmp)+1],22*[cos(tmp),sin(tmp)]);

axis off equal
set(gcf, 'menu', 'none', 'NumberTitle','off','name','EMDLAB')
p = patch('faces',m.m.cl, 'vertices', m.m.nodes, 'faceColor', 'c');
for i = 1:1000
m.rotateInner(pi/30);
m.m.Ne
p.Faces = m.m.cl;
p.Vertices = m.m.nodes;
drawnow
pause(0.02);
end
