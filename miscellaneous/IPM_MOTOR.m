%% initialization
clc
clear
materialdir = [cd,'\MaterialsData'];
g = GDBC2D;

p.h = 4;
p.wy = 1.5;
p.wc = 1;
p.g = 0.2;

g = g.newdlinewdkps([-p.h/2-p.wy,0],[-p.h/2,0],'Nnodes',5);
g = g.newdlinewdkps([-p.h/2,0],[-p.h/2,p.h/2],'Nnodes',6);
g = g.newdlinewdkps([-p.h/2,p.h/2],[-p.h/2+p.wc,p.h/2],'Nnodes',4);
g = g.newdlinewdkps([-p.h/2+p.wc,p.h/2],[p.h/2,p.h/2],'Nnodes',8);
g = g.newdlinewdkps([p.h/2,p.h/2],[p.h/2,p.g/2],'Nnodes',6);
g = g.newdlinewdkps([p.h/2,p.g/2],[p.h/2+p.wy,p.g/2],'Nnodes',5);
g = g.newdlinewdkps([p.h/2+p.wy,p.g/2],[p.h/2+p.wy,p.h/2+p.wy],'Nnodes',9);
g = g.newdlinewdkps([p.h/2+p.wy,p.h/2+p.wy],[-p.h/2-p.wy,p.h/2+p.wy],'Nnodes',14);
g = g.newdlinewdkps([-p.h/2-p.wy,p.h/2+p.wy],[-p.h/2-p.wy,p.h/2],'Nnodes',5);
g = g.newdlinewdkps([-p.h/2-p.wy,p.h/2],[-p.h/2-p.wy,0],'Nnodes',6);

g = g.newdlinewdkps([-p.h/2,0],[-p.h/2+p.wc,0],'Nnodes',4);
g = g.newdlinewdkps([-p.h/2+p.wc,0],[-p.h/2+p.wc,p.h/2],'Nnodes',6);
g = g.cmirrorbd('L3',[-p.h/2-p.wy/2,0],[-p.h/2-p.wy/2,1]);
g = g.cmirrorbd('L12',[-p.h/2-p.wy/2,0],[-p.h/2-p.wy/2,1]);
g = g.cmirrorbd('L11',[-p.h/2-p.wy/2,0],[-p.h/2-p.wy/2,1]);

g = g.newdlinewdkps([-p.h/2+p.wc,0],[0,0],'Nnodes',4);
g = g.newdlinewdkps([0,0],[p.h/2+p.wy+p.wc,0],'Nnodes',12);
g = g.newdlinewdkps([p.h/2+p.wy+p.wc,0],[3*p.h,0],'Nnodes',10,'alpha',1.2);
g = g.newdlinewdkps([3*p.h,0],[3*p.h,3*p.h],'Nnodes',8);
g = g.newdlinewdkps([3*p.h,3*p.h],[-3*p.h,3*p.h],'Nnodes',15);
g = g.newdlinewdkps([-3*p.h,3*p.h],[-3*p.h,0],'Nnodes',8);
g = g.newdlinewdkps([-p.h/2-p.wy-p.wc,0],[-3*p.h,0],'Nnodes',10,'alpha',1.2);

g = g.newcb('y','L1',1,'L2',1,'L3',1,'L4',1,'L5',1,'L6',1 ...
    ,'L7',1,'L8',1,'L9',1,'L10',1);
g = g.newcb('c1','L11',1,'L12',1,'L3',-1,'L2',-1);
g = g.newcb('c2','L13',1,'L14',-1,'L15',-1,'L10',-1);
g = g.newcb('a','L16',1,'L17',1,'L18',1,'L19',1 ...
    ,'L20',1,'L21',1,'L22',-1,'L14',1 ...
    ,'L13',-1,'L9',-1,'L8',-1,'L7',-1 ...
    ,'L6',-1,'L5',-1,'L4',-1,'L12',-1 );

g = g.newdDM('y','y');
g = g.newdDM('c1','c1');
g = g.newdDM('c2','c2');
g = g.newdDM('a','a');

m = MDBCT(g);clear g;

m = m.setMaterial('y','m19_24ga');
m = m.setMaterial('c1','copper');
m = m.setMaterial('c2','copper');
m = m.addMaterial(materialdir,'air');
m = m.addMaterial(materialdir,'m19_24ga');
m = m.addMaterial(materialdir,'copper');

m = m.cmirrormz('y2','y',[1,0]);
m = m.cmirrormz('a2','a',[1,0]);
m = m.cmirrormz('c22','c2',[1,0]);
m = m.cmirrormz('c11','c1',[1,0]);

m = m.joinmzs('yoke','y','y2');
m = m.joinmzs('air','a','a2');
m = m.joinmzs('coil1','c1','c11');
m = m.joinmzs('coil2','c2','c22');

m = m.ggmesh;
m.showmzs
%% solver setting
s = IHNLNRMSTL3(m);clear m;
% setting units
s.scs.l = 1e-3;
s.scs.f = 1e6;
%% proccess

s = s.setExcitation('coil1',-10);
s = s.setExcitation('coil2',10);

s.m = s.m.ggmesh;
s.m = s.m.strefine;
s.m = s.m.strefine;
k0 = s.m.getfb;

s.m = s.m.evalKeFeC('TL3');
s = s.clearallbcs;
s = s.setdbc(unique(k0(:)),0);
s = s.assignEdata;
s = s.solve(1e-8,30);

s.plotBmagw('yoke')

s.evalLF('coil1')-s.evalLF('coil2')