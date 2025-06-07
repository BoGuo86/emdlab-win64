
clc
clear

global gdb mdb sdb
gdb = G_gdbc();
mdb = mdbc();
sdb = sdbc();

h1 = 0.25;
h2 = 0.25;

newcirculard('cr1',[1,-1],0.25,0.9*h1)
newcirculard('cr2',[1,0],0.25,0.9*h1)
newcirculard('cr3',[1,1],0.25,0.9*h1)

newcirculard('cl1',[-1,-1],0.25,0.9*h1)
newcirculard('cl2',[-1,0],0.25,0.9*h1)
newcirculard('cl3',[-1,1],0.25,0.9*h1)

newcircularcb('c',[0,0],10,0.9*h2)

newd('a','c',{'CB1','CB2','CB3','CB4','CB5','CB6'})

[p,t,ss,ipos] = gmesh('cr1',h1);
newmz('cr1m',p,t,ss,ipos);

[p,t,ss,ipos] = gmesh('cr2',h1);
newmz('cr2m',p,t,ss,ipos);

[p,t,ss,ipos] = gmesh('cr3',h1);
newmz('cr3m',p,t,ss,ipos);

[p,t,ss,ipos] = gmesh('cl1',h1);
newmz('cl1m',p,t,ss,ipos);

[p,t,ss,ipos] = gmesh('cl2',h1);
newmz('cl2m',p,t,ss,ipos);

[p,t,ss,ipos] = gmesh('cl3',h1);
newmz('cl3m',p,t,ss,ipos);

[p,t,ss,ipos] = gmesh('a',h2);
newmz('am',p,t,ss,ipos);

ggmesh

tic
am2mdb('air')
setff('cr1m',-1)
setff('cr2m',-1)
setff('cr3m',-1)
setff('cl1m',1)
setff('cl2m',1)
setff('cl3m',1)
FMSEqs()
setdbcontotal(0)
solveKUF
toc

evalgradu()
plotmag