
f = [1,2,3;2,4,3];
p = [0,0;1,0;0,1;1,1];
v = [1;0;0;-3];

c = tmzpc_contour_tl3(f,p,v,linspace(min(v),max(v),50));

hold all
patch('faces', f, 'vertices', p,'facecolor', 'w');

t = 1:size(c,1);
t = reshape(t,2,[])';
patch('faces', t, 'vertices', c,'facecolor', 'b');

axis off equal
