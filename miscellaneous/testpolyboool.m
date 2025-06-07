
clc
clear all

global gdb mdb

gdb = gdbc;
mdb = mdbc;
h = 0.1;
newrectangulard('r1',0,0,1,1,0.9*h,0.9*h)
newrectangulard('r2',0.2,0.2,0.6,0.6,0.9*h,0.9*h)
dsubtract1('a','r1','r2')

newrectangulard('r3',-0.5,-0.5,1,1,0.9*h,0.9*h)
dunion1('b','a','r3')


[p,t,ss,ipos] = gmesh('b',h);
newmz('am',p,t,ss,ipos);

ggmesh
plotgm
