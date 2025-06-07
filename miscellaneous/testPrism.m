p = [0,1;0,0;1,0];

p = [p,[0;0;0];p,[1;1;1]];

t = [6,4,3,5
    5,3,2,1
    5,4,3,1];

subplot(121)
tr= triangulation(t,p);
tetramesh(tr)

subplot(122)
patch('Faces',tr.freeBoundary,'Vertices',...
tr.Points,'FaceColor','c','faceAlpha',0.5)
axis equal off
rotate3d
