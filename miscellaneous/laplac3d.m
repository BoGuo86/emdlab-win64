% k = 50;
% p = 2*rand(k,3)-1;

[x,y,z] = meshgrid(-1:0.2:1);
p = unique([x(:),y(:),z(:)],'rows');

index = sqrt(sum(p.^2,2))<0.9;
p = p(index,:);
k = size(p,1);


[x,y,z] = sphere(20);
p1 = unique([x(:),y(:),z(:)],'rows');
% k = size(p,1);

% p1 = [0,0,0
%     0,0,1
%     0,1,0
%     0,1,1
%     1,0,0
%     1,0,1
%     1,1,0
%     1,1,1];

% p1 = [p1' ; (p1+circshift(p1,-1))'/2];
% p1 = reshape(p1,2,[]);
% p1 = p1';
% 
% p1 = [p1' ; (p1+circshift(p1,-1))'/2];
% p1 = reshape(p1,2,[]);
% p1 = p1';

p = [p;p1];

t = delaunayn(p);
tr = triangulation(t,p);

% subplot(121)
%  tetramesh(tr,'facecolor','none')
% axis equal off
% title('Initial Mesh')

c = tr.edges;
n = max(max(c));
c = sparse(c(:,1),c(:,2),ones(size(c,1),1),n,n);
amp = 1./sum(c(1:k,:),2);
amp = sparse(1:length(amp),1:length(amp),amp);
c = c + c';

clc

% subplot(122)
for i = 1:100
pnew = c(1:k,:)*p;
pnew = diag(1./sum(c(1:k,:),2))*pnew;

d = sqrt(sum((p(1:k,:)-pnew).^2,2));
disp(sum(d));
p(1:k,:) = pnew;

% tr = triangulation(t,p);
%  tetramesh(tr,'facecolor','none')

% axis equal off
% pause(0.01)

if d < 0.001
    disp(i)
    break
end
    
end
title('Final Mesh')

c = tr.incenter;

close all
index = bitand(c(:,1)>0,c(:,2)>0);
index = bitand(index,c(:,3)>0);
t = t(index,:);
tr = triangulation(t,p);
face = tr.freeBoundary;
patch('faces',face,'vertices',p,'facecolor','c')
 axis off equal
 view([-1,-1,-1])
