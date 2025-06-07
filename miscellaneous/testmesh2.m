clc
clear

h = 0.1;
p1 = [0,0];
p2 = [2,0];
p3 = [2,1];
p4 = [1,1];
p5 = [1,2];
p6 = [0,2];

p = [p1;
    getlinepoints(p1,p2,h)
    p2;
    getlinepoints(p2,p3,h)
    p3;
    getlinepoints(p3,p4,h)
    p4
    getlinepoints(p4,p5,h)
    p5
    getlinepoints(p5,p6,h)
    p6
    getlinepoints(p6,p1,h)];

[x,y] = meshgrid(min(p(:,1)):h:max(p(:,1)),min(p(:,2)):h*sqrt(3)/2:max(p(:,2)));
x(2:2:end,:)=x(2:2:end,:)+h/2;  
ip = [x(:),y(:)];

%  ip = rand(20,2);

[s,t] = inpolygon(ip(:,1),ip(:,2),p(:,1),p(:,2));
ip = ip(s&~t,:);

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
    if mov<0.01
        disp(iter)
        break
    end
end


c = dt.incenter;
s = inpolygon(c(:,1),c(:,2),p(:,1),p(:,2));
t = dt.ConnectivityList(s,:);

triplot(triangulation(t,[p;ip]));axis off equal

