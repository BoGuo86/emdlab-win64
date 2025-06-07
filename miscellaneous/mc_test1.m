
x = linspace(0,5,12);

p1 = [x;1*x]';
p2 = [x;1*x+0.5]';

m = mc_lineAirGap(p1,p2);


axis off equal
p = patch('faces',m.m.cl, 'vertices', m.m.nodes, 'faceColor', 'c');
for i = 1:10
m.shiftMovingPoints(0.2);
m.m.Ne
p.Faces = m.m.cl;
p.Vertices = m.m.nodes;
drawnow
pause(0.5);
end
