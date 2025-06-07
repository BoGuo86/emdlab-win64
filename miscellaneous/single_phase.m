
tic
clc
clear all

global gdb mdb sdb

gdb = G_gdbc;
mdb = M_mdbc;
sdb = MSTL1_sdbc;

h = 0.25;
hy = 0.25;
a = 2;
b = 2;
wy = 1;
wc = 0.5;
g = 0.3;

newpolygonald('yoke',[-a-wy,0;-a,0;-a,b;-a+wc,b;a,b;a,g/2;...
    a+wy,g/2;a+wy,b+wy;-a-wy,b+wy;-a-wy,b],0.9*hy)
newdlinewdkps([-a,0],[-a+wc,0],0.9*h)
newdlinewdkps([-a+wc,0],[-a+wc,b],0.9*h)
newdlinewdkps([-a-wy,0],[-a-wy-wc,0],0.9*h)
newdlinewdkps([-a-wy-wc,0],[-a-wy-wc,b],0.9*h)
newdlinewdkps([-a-wy-wc,b],[-a-wy,b],0.9*h)

newcbwd('c1',{'L11','L12','L3','L2'},[1 1 -1 -1])
newcbwd('c2',{'L13','L14','L15','L10'},[1 1 1 1])

scale = 5;
newdlinewdkps([-a+wc,0],[0,0],0.9*h)
newdlinewdkps([0,0],[scale*a,0],0.9*h)
newdlinewdkps([scale*a,0],[scale*a,scale*b],0.9*h)
newdlinewdkps([scale*a,scale*b],[-scale*a,scale*b],0.9*h)
newdlinewdkps([-scale*a,scale*b],[-scale*a,0],0.9*h)
newdlinewdkps([-scale*a,0],[-a-wy-wc,0],0.9*h)

newcbwd('air',{'L16','L17','L18','L19',...
    'L20','L21','L14','L15',...
    'L9','L8','L7','L6',...
    'L5','L4','L12'
    },[1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1])


[p,t,ss,ipos] = gmesh('c1',hy);
newmz('c1m',p,t,ss,ipos);

[p,t,ss,ipos] = gmesh('c2',hy);
newmz('c2m',p,t,ss,ipos);

[p,t,ss,ipos] = gmeshMG1('air',h);
newmz('airm',p,t,ss,ipos);

[p,t,ss,ipos] = gmesh('yoke',hy);
newmz('ym',p,t,ss,ipos,'m19_24ga');

ggmeshnew
plotgm
%%
am2mdb('air')
am2mdb('iron')
am2mdb('m19_24ga')
am2mdb('copper')
evalKeFeCnew()
setddbc('airm_L18',0)
setddbc('airm_L19',0)
setddbc('airm_L20',0)
%%
% c = 0:10:50;
% 
% y = c;
% for i = 2:length(c)
setff('c1m',-2000)
setff('c2m',2000)


% setting scalers
sdb.scs.l = 1e-3;
sdb.scs.k = sdb.pcts.nu0;
sdb.scs.u = 1/sdb.pcts.nu0;
sdb.scs.f = 1e6;
IHNL_NR_MS_TL1_solver(1e-3,1e-2,20,0.01)
% LsolverMSd1()
% LsolverMSd1
% IHL_MS_TL1_solver
% evalgradu()
% evalEnorm
% plotmag
% plotucontour()
% b = sqrt(sum(sdb.B.^2,2));
% y(i) = -ms_evalLF('c1m')+ms_evalLF('c2m');
% 
% end
% ms_evalLF('c1m')+ms_evalLF('c2m')
% % ms_evalB
ms_plotBmagSmooth
% 
% ms_plotBmag
% close
% plot(c,y)
