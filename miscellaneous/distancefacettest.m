
p = rand(200,3);
p1 = [1,0,0];
p2 = [0,1,0];
p3 = [0,0,1];


trisurf([1,2,3],p1,p2,p3,'facecolor','y');hold on

d = sdpsffacet(p,p1,p2,p3);
plot3(p(d>0,1),p(d>0,2),p(d>0,3),'.','color','b');axis off equal;hold on
plot3(p(d<0,1),p(d<0,2),p(d<0,3),'.','color','r');hold on
plot3(p(d==0,1),p(d==0,2),p(d==0,3),'.','color','g');