k = 50;
p = rand(k,2);

% [x,y] = meshgrid(0.2:0.1:0.8);
% p = unique([x(:),y(:)],'rows');
% k = size(p,1);

p1 = [0,0;1,0;1,1;0,1];
p1 = [p1' ; (p1+circshift(p1,-1))'/2];
p1 = reshape(p1,2,[]);
p1 = p1';

p1 = [p1' ; (p1+circshift(p1,-1))'/2];
p1 = reshape(p1,2,[]);
p1 = p1';
p = [p;p1];

t = delaunay(p);
tr = triangulation(t,p);

subplot(121)
triplot(tr)
axis equal off
title('Initial Mesh')

c = tr.edges;
n = max(max(c));
c = sparse(c(:,1),c(:,2),ones(size(c,1),1),n,n);
c = c + c';

clc

subplot(122)
for i = 1:100
pnew = c(1:k,:)*p;
pnew = diag(1./sum(c(1:k,:),2))*pnew;

d = sqrt(sum((p(1:k,:)-pnew).^2,2));
disp(sum(d));
p(1:k,:) = pnew;

tr = triangulation(t,p);
triplot(tr)

axis equal off
pause(0.01)

if d < 0.001
    disp(i)
    break
end
    
end
title('Final Mesh')
