function mz = STLBIN2TTMZ(Dir)

f = fopen('C:\Users\AliJamalifard\Desktop\Part1.STL','r');

fread(f,80,'*char');

Nt = fread(f,1,'uint32');

n = zeros(3*Nt,1);
p = zeros(9*Nt,1);

for i = 1:Nt
    n(3*i-2:3*i) = fread(f,3,'single');
    p(9*i-8:9*i-6) = fread(f,3,'single');
    p(9*i-5:9*i-3) = fread(f,3,'single');
    p(9*i-2:9*i) = fread(f,3,'single');
    fread(f,1,'uint16');
end

fclose(f);

p = reshape(p,3,[]);
p = p';

t = 1:3*Nt;
t = reshape(t,3,[]);
t = t';

[p,~,index] = uniquetol(p,1e-4,'ByRows',true);

t = index(t);

fd = @(x) dpoly3d(x,t,p);

mz = G3DFV(t,p);
% h0 = 1.5;
% [p,t] = distmeshnd(fd,@huniform,1.05*h0,[min(p);max(p)],[p;getbn(G3DFV(t,p),h0)]);
%  
% % t = delaunayn(p);
% % c = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:)+p(t(:,4),:))/4;
% % d = fd(c);
% % t = t(d<-1e-6,:);
% mz = TTMDBC({'a',TTMZPC(t,p)});
% close all
% mz.showmzs


end