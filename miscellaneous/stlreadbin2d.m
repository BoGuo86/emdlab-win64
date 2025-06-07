
f = fopen('C:\Users\AliJamalifard\Desktop\b.STL','r');

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

hold on;
trisurf(index(t),p(:,1),p(:,2),p(:,3),'FaceColor','c','edgeColor','w');

set(gcf,'Color','k')
 set(gcf,'Renderer','opengl');
axis equal off

view([0 0 1])
%%
% tr = triangulation(index(t),p(:,1:2));
% 
% triplot(tr);axis off equal

      