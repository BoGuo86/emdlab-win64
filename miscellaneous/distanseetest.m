
p = rand(5000,2);

tic
d = sdfpolygon(p,[0.25,0.25;0.75,0.25;0.25,0.75]);
toc
tic
d = dpoly(p,[0.25,0.25;0.75,0.25;0.25,0.75]);
toc
plot(p(d>0,1),p(d>0,2),'.','color','r');hold on
plot(p(d<0,1),p(d<0,2),'.','color','b')