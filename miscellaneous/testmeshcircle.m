clc
clear

h = 0.1;
p1 = [1,0];
p2 = [-1,0];

p = [p1;
    getarcpointscpp([0,0],p1,p2,1,h)
    p2;
    getarcpointscpp([0,0],p2,p1,1,h)];

[x,y] = meshgrid(min(p(:,1)):h:max(p(:,1)),min(p(:,2)):h*sqrt(3)/2:max(p(:,2)));
x(2:2:end,:)=x(2:2:end,:)+h/2;  
ip = [x(:),y(:)];

[s,t] = inpolygon(ip(:,1),ip(:,2),p(:,1),p(:,2));
ip = ip(s&~t,:);

% l = sqrt(sum(p.^2,2));
% ip = [diag(0.8*l)*p;diag(0.6*l)*p;diag(0.4*l)*p;
%     diag(0.2*l)*p;[0,0]];


dt = delaunayTriangulation([p;ip]);
triplot(dt);axis off equal

Np = size(p,1);
Nip = size(ip,1);

t = cell(1,size(ip,1));
index = 1:Nip;

for iter = 1:100
    mov = 0;
    v = dt.vertexAttachments((1+Np:Np+Nip)');
    for i = 1:size(ip,1)
        t{i} = dt.ConnectivityList(v{i},:);
        t{i} = t{i}(:);
        t{i} = unique(t{i});
        t{i} = t{i}(t{i}~=i+Np);
    end
    for i = index
        pp = dt.Points(t{i},:);
        pp = sum(pp)/length(t{i});
        movi = norm(dt.Points(i+Np,:)-pp);
        mov = mov + movi;
        if movi<0.001
            index = setdiff(index,i);
        end
        ip(i,:) = pp;
        dt.Points(Np+i,:) = pp;
    end
%     dt.Points(1+Np:Np+Nip,:) = ip;
    cla,triplot(dt);axis off equal
    drawnow
    if mov<0.001
        disp(iter)
        break
    end
end

