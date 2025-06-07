
clc
clear

c1 = get_circle(5);
c2 = get_circle(7);
c3 = get_circle(3.5);

n = 6;
r(1:n) = polyshape([3,6,6,3],[-1,-1,1,1]);
for i = 1:n
    r(i).Vertices = ext_protate2(r(i).Vertices, (i-1)*2*pi/n);
end

s = subtract(c2,c1);


for i = 1:n
    s = union(s,r(i));
end
s = subtract(s,c3);

c = get_circle(10);
c = subtract(c,s);

m = TMDBC;
t = triangulation(s);
m.addmz('s', TMZPC(t.ConnectivityList, t.Points));

t = triangulation(c);
m.addmz('c', TMZPC(t.ConnectivityList, t.Points));
m.showg

function p = get_circle(r)
t = linspace(0,2*pi,50);
t = t(1:end-1);
x = r*cos(t);
y = r*sin(t);
p = polyshape(x,y);
end
