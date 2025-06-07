tic
clc
clear

global gdb mdb sdb

gdb = gdbc;
mdb = mdbc;
sdb = sdbc;

% main dim
rsh = 8.000000;
R0 = 30.400000;
R1 = 40.900000;
R2 = 65.800000;
R3 = 74.900000;
betar = 19 * pi/180;
betas = 16 * pi/180;
g = 0.3; 
dw = 4;
Rw = R1+g+dw;

gammar = asin(R1*sin(betar/2)/R0);
gammas = asin((R1+g)*sin(betas/2)/R2);
gammaw = asin((R1+g)*sin(betas/2)/Rw);

p1 = R0*[cos(gammar),sin(-gammar)];
p2 = R1*[cos(betar/2),sin(-betar/2)];
p3 = R1*[cos(betar/2),sin(betar/2)];
p4 = R0*[cos(gammar),sin(gammar)];


newkp('p1',[rsh,0])
newkp('p2',[R1,0])
newkp('p3',p3)
newkp('p4',p4)
newkp('p5',R0*[cos(pi/8),sin(pi/8)])
newkp('p6',rsh*[cos(pi/8),sin(pi/8)])

hr = 2;
hs =2;
hg =1;
thetag =0.5;

newline('l1','p1','p2',0.9*hr)
newarc('l2',[0,0],'p2','p3',0.9*min(hr,hg),thetag)
newline('l3','p3','p4',0.9*min(hr,hg))
newarc('l4',[0,0],'p4','p5',0.9*0.9*min(hr,hg),5)
newline('l5','p6','p5',0.9*hr)
newarc('l6',[0,0],'p1','p6',0.9*hr,5)




newcb('srb',getlist('l',1:6),[1 1 1 1 -1 -1])
newd('srotor','srb')



sp = [R1+g,0;R3,0;R3*[cos(pi/12),sin(pi/12)];...
    R2*[cos(pi/12),sin(pi/12)];R2*[cos(gammas),sin(gammas)]
    Rw*[cos(gammaw),sin(gammaw)]
    (R1+g)*[cos(betas/2),sin(betas/2)]];
newdkps(sp);



newline('l7','kp1','kp2',0.9*hs)
newarc('l8',[0,0],'kp2','kp3',0.9*hs,10)
newline('l9','kp4','kp3',0.9*hs)
newarc('l10',[0,0],'kp5','kp4',0.9*min(hs,hg),5)
newline('l11','kp5','kp6',0.9*min(hs,hg))
newline('l12','kp6','kp7',0.9*min(hs,hg))
newarc('l13',[0,0],'kp1','kp7',0.9*min(hs,hg),thetag)


newcb('ssb',getlist('l',7:13),[1 1 -1 -1 1 1 -1])
newd('sstator','ssb')



newkp('wkp',Rw*[cos(pi/12),sin(pi/12)])
newline('L1','wkp','kp4',0.9*hg)
newarc('L2',[0,0],'kp6','wkp',0.9*hg,thetag)

newcb('swb',{'L2','l11','l10','L1'},[-1 -1 1 -1])
newd('sw','swb')



newdkp([(R1+g/2),0])
newdkp((R1+g/2)*[cos(pi/4),sin(pi/4)])

cmirrors('a1','l2',[cos(pi/8),sin(pi/8)])
cmirrors('a2','l3',[cos(pi/8),sin(pi/8)])
cmirrors('a3','l4',[cos(pi/8),sin(pi/8)])
newarc('ra',[0,0],'p3','kp11',0.9*hg,thetag)

newcb('rab',{'l3','l4','a3','a2','ra'},[1 1 -1 -1 -1])
newd('rap1','rab')


cmirrors('a4','l13',[cos(pi/12),sin(pi/12)])
cmirrors('a5','l12',[cos(pi/12),sin(pi/12)])
cmirrors('a6','L2',[cos(pi/12),sin(pi/12)])

newkp('st',(R1+g)*[cos(pi/12),sin(pi/12)])
newarc('sa1',[0,0],'kp7','st',0.9*hg,thetag)
% newarc('sa2',[0,0],'st','kp14',0.9*hg,thetag)
cmirrors('sa2','sa1',[cos(pi/12),sin(pi/12)])


cmirrors('x1','a4',[cos(pi/6),sin(pi/6)])
cmirrors('x2','sa2',[cos(pi/6),sin(pi/6)])

newcb('sab',{'l12','L2','a6','a5','sa2','sa1'},[-1 1 -1 1 1 -1])
newd('sap1','sab')

newarc('ga',[0,0],'kp8','kp9',0.9*hg,thetag)

newline('gs1','p2','kp8',0.9*hg)
newline('gs2','kp8','kp1',0.9*hg)
newline('gs3','kp10','kp9',0.9*hg)
newline('gs4','kp9','kp17',0.9*hg)

newcbwd('rag',{'gs1','ga','gs3','a1','ra','l2'},[1 1 -1 1 -1 -1])

newcbwd('sag',{'gs2','l13','sa1','sa2','a4','x1','x2','gs4','ga'},[1 1 1 -1 -1 1 1 -1 -1])

%%
% meshing
% rotor

[p,t,ss,ipos] = gmesh('srotor',hr);
newmz('rm1',p,t,ss,ipos,'iron');
cmirrormz('rm2','rm1',[cos(pi/8),sin(pi/8)])
crotatemz('rm3','rm1',pi/4)
crotatemz('rm4','rm2',pi/4)

% rotor air pocket
[p,t,ss,ipos] = gmesh('rap1',hg);
newmz('rapm1',p,t,ss,ipos,'air');
crotatemz('rapm2','rapm1',pi/4)

% stator mesh
[p,t,ss,ipos] = gmesh('sstator',hs);
newmz('sm1',p,t,ss,ipos,'iron');
cmirrormz('sm2','sm1',[cos(pi/12),sin(pi/12)])
crotatemz('sm3','sm1',pi/6)
crotatemz('sm4','sm2',pi/6)
crotatemz('sm5','sm3',pi/6)
crotatemz('sm6','sm4',pi/6)

% stator air pocket
[p,t,ss,ipos] = gmesh('sap1',hg);
newmz('sapm1',p,t,ss,ipos,'air');
crotatemz('sapm2','sapm1',pi/6)
crotatemz('sapm3','sapm2',pi/6)



[p,t,ss,ipos] = gmesh('sw',hg);
newmz('c1a1m',p,t,ss,ipos,'air');
cmirrormz('c2a2m','c1a1m',[cos(pi/12),sin(pi/12)])
crotatemz('c2a1m','c1a1m',pi/6)
crotatemz('c3a2m','c2a2m',pi/6)
crotatemz('c3a1m','c2a1m',pi/6)
crotatemz('c1a2m','c3a2m',pi/6)


% newdlinewdkps([0,0],[rsh,0],0.9*hr)
% newdlinewdkps([0,0],rsh*[cos(pi/8),sin(pi/8)],0.9*hr)
% newcb('shb',{'L3','l6','L4'},[1 -1 -1])
% newd('shaft','shb')
% 
% [p,t,ss,ipos] = gmesh('shaft',hr);
% newmz('shm',p,t,ss,ipos);
% 
% newdlinewdkps([R3,0],[R3+4,0],0.9*hs)
% newdlinewdkps(R3*[cos(pi/12),sin(pi/12)],(R3+4)*[cos(pi/12),sin(pi/12)],0.9*hs)
% 
% newarc('ha',[0,0],'kp19','kp20',5,0.9*hs)
% 
% newcb('hb',{'L5','ha','L6','l8'},[1 1 -1 -1])
% newd('h','hb')
% 
% [p,t,ss,ipos] = gmesh('h',hs);
% newmz('hm',p,t,ss,ipos);

[p,t,ss,ipos] = gmesh('rag',hg);
newmz('rag1m',p,t,ss,ipos);
cmirrormz('rag2m','rag1m',[cos(pi/4),sin(pi/4)])

[p,t,ss,ipos] = gmesh('sag',hg);
newmz('sag1m',p,t,ss,ipos);
cmirrormz('sag2m','sag1m',[cos(pi/4),sin(pi/4)])

ggmesh

toc

%%
tic
am2mdb('air')
am2mdb('iron')
am2mdb('copper')
% setopbc('rm1_l1', 'rm4_l1', 1)
% setopbc('sm1_l7', 'sm6_l7', 1)
% setopbc('rag1m_gs1', 'rag2m_gs1',1)
% setopbc('sag1m_gs2', 'sag2m_gs2', 1)
% 
% setdbc('rm1_l6',0)
% setdbc('rm2_l6',0)
% setdbc('rm3_l6',0)
% setdbc('rm4_l6',0)
% 
% setdbc('sm1_l8',0)
% setdbc('sm2_l8',0)
% setdbc('sm3_l8',0)
% setdbc('sm4_l8',0)
% setdbc('sm5_l8',0)
% setdbc('sm6_l8',0)
evalsijfid1()
setff('c1a1m',0.01)
setff('c1a2m',0.01)
LsolverMSd1()
setdbcontotal(0)
solveKUF
toc

evalgradu()
evalEnorm
% plotmag
plotucontour()



