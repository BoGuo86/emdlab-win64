function y = pole_shape(y)


materialdir = [cd,'\MaterialsData'];
g = GDBC2D;

p.g = 0.5;
p.Wp = 10;
p.Wpsh = 8;
p.Wmid = 8;
p.Hp = 30;
p.Hc = 20;

% y = [0,y,1];
y = y';
Ny = length(y);
tmp = linspace(0,p.Wp+p.Wpsh,Ny);



tmp = [tmp',p.Hp-y];
for i = 1:Ny-1
    g = g.newdlinewdkps(tmp(i,:),tmp(i+1,:),'maxLength',1);
end
g = g.newdlinewdkps(tmp(end,:),[p.Wp+p.Wpsh,p.Hc],'maxLength',1);
g = g.newdlinewdkps([p.Wp+p.Wpsh,p.Hc],[p.Wp,p.Hc],'maxLength',1);
g = g.newdlinewdkps([p.Wp,p.Hc],[p.Wp,0],'maxLength',1);
g = g.newdlinewdkps([p.Wp,0],[p.Wp+p.Wpsh,0],'maxLength',1);
g = g.newdlinewdkps([p.Wp+p.Wpsh,0],[p.Wp+p.Wpsh+p.Wmid,0],'maxLength',1);
g = g.newdlinewdkps([p.Wp+p.Wpsh+p.Wmid,0],[p.Wp+p.Wpsh+p.Wmid,p.Hp+p.g],'maxLength',1);
g = g.newdlinewdkps([p.Wp+p.Wpsh+p.Wmid,p.Hp+p.g],[0,p.Hp+p.g],'maxLength',1);
g = g.newdlinewdkps([0,p.Hp+p.g],tmp(1,:),'maxLength',1);
g = g.newdlinewdkps([p.Wp+p.Wpsh,0],[p.Wp+p.Wpsh,p.Hc],'maxLength',1);
tmp1 = [p.Wp+p.Wpsh+p.Wmid,p.Hp+p.g];
tmp2 = [0,p.Hp+p.g];
g = g.newdlinewdkps(tmp1,tmp1+[0,p.g],'maxLength',1);
g = g.newdlinewdkps(tmp1+[0,p.g],tmp2+[0,p.g],'maxLength',1);
g = g.newdlinewdkps(tmp2+[0,p.g],tmp2,'maxLength',1);

a = cell(1,2*Ny);
for i = 1:Ny
    a{2*i-1} = ['L',num2str(i)];
    a{2*i} = 1;
end
g = g.newcb('a',a{:} ...
    ,['L',num2str(Ny+8)],-1 ...
    ,['L',num2str(Ny+4)],1 ...
    ,['L',num2str(Ny+5)],1 ...
    ,['L',num2str(Ny+6)],1 ...
    ,['L',num2str(Ny+7)],1);

g = g.newcb('c' ...
    ,['L',num2str(Ny+8)],1 ...
    ,['L',num2str(Ny+1)],1 ...
    ,['L',num2str(Ny+2)],1 ...
    ,['L',num2str(Ny+3)],1);

g = g.newcb('g' ...
    ,['L',num2str(Ny+9)],1 ...
    ,['L',num2str(Ny+10)],1 ...
    ,['L',num2str(Ny+11)],1 ...
    ,['L',num2str(Ny+6)],-1);

g = g.newdd('a',1,'a');
g = g.newdd('c',1,'c');
g = g.newdd('g',1,'g');

m = MDBCT(g);

m = m.addMaterial(materialdir,'air');

% m = m.cmirrormz('g1','g',[0,1]);
% m = m.cmirrormz('c1','c',[0,1]);
% m = m.cmirrormz('a1','a',[0,1]);
%% solver setting
s = IHNLNRMSTL3(m);clear m;
% setting units
s.scs.l = 1e-3;
s.scs.f = 1e6;
%% proccess

s = s.setExcitation('c',-1);
% % s = s.setExcitation('c1',10);
% 
s.m = s.m.ggmesh;
s.m = s.m.strefine;
% s.m = s.m.strefine;
% % s.m = s.m.strefine;
% 
% % k = [s.m.getIndexOnRay([-p.W,0],[-p.W,1])
% %     s.m.getIndexOnRay([p.W,0],[p.W,1])];
% 
k0 = s.m.getIndexOnRay([0,0],[0,1]);
% 
% % [km,ks] = s.m.splitShift(k,[2*p.W,0]);
% s.m.showmeshfb;
% hold on;plot(s.m.p(k0,1),s.m.p(k0,2),'*','color','g');
% % hold on;plot(s.m.p(km,1),s.m.p(km,2),'*','color','g');
% % hold on;plot(s.m.p(ks,1),s.m.p(ks,2),'*','color','g');
% 
s.m = s.m.evalKeFeC('TL3');
s = s.clearallbcs;
s = s.setdbc(unique(k0(:)),0);
% s = s.setopbc(km,ks);
s = s.assignEdata;
s = s.solve(1e-8,30);


y = s.getbonb('g',['L',num2str(Ny+6)]);
y = y(:,2);
% y = y';
% x = (0:30);
% 
% a = [1,3,5,7,9,11,13,15];
% Na = length(a);
% a = repmat(a',1,length(x));
% 
% c = cos(a.*repmat(x,Na,1));
% c = c.*(repmat(y,Na,1));
% c = sum(c,2);

% y = sqrt(sum(c(2:end).^2))/c(1);

y = flipud(y);
y = [y(1:end-1);-flipud(y)];
y = [y(1:end-1);flipud(y)];
% plot(y)
y = fft(y);
y = 2*real(y(2:(length(y)-1)/2))/length(a);

clc
y = sqrt(sum(y(2:end).^2))/y(1)

close all
s.m.plotwf
pause(0.1)

end



