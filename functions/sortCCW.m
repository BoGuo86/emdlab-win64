function [p,index] = sortCCW(p,ali)

if nargin<2
    ali = 1;
end

Np = size(p,1);
index = 1:Np;
index = index';

for i = 1:Np-1
    for j = 1:Np-i
        if p(j,1)*p(j+1,2)-p(j,2)*p(j+1,1)<0
            temp = p(j,:);
            p(j,:) = p(j+1,:);
            p(j+1,:) = temp;
            temp = index(j);
            index(j) = index(j+1);
            index(j+1) = temp;
        end
    end
end

if ali
    while det([p(1,:);p(end,:)])<-1e-3
        %     det([p(1,:);p(end,:)])
        %     plot(p(:,1),p(:,2))
        %     axis off
        %     pause(0.25)
        p = circshift(p,1);
        index = circshift(index,1);
    end
end

end