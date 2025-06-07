x = linspace(-1,1,1000);
% y = @(x) sin(pi*x)-0.5*sin(5*pi*x)+0.05*sin(7*pi*x);
y = @(x) exp(-100*abs(x));
xi = linspace(-1,1,5);
yi = feval(y,xi);

n = length(xi);
Intervals = [1:n-1;2:n]';
Flags = ones(1,n-1);

plot(x,feval(y,x));hold on
plot(xi,yi);hold off
axis off

for iter = 1:100
    
    n = size(Intervals,1);
    if ~Flags
%         clc
        disp(iter)
        break
    end
        
    
    for i = 1:n
        
        if Flags(i)
            
            xtmp = mean(xi(Intervals(i,:)));
            ytmp = mean(yi(Intervals(i,:)));
            
            ymid = feval(y,xtmp);
            
            l1 = getL1(xi(Intervals(i,1)),0,xi(Intervals(i,2)),ymid-ytmp);
            l2 = getL1(xi(Intervals(i,1)),ymid-ytmp,xi(Intervals(i,2)),0);
            
            L = abs(ymid-ytmp)*diff(xi(Intervals(i,:)))
            
            if L>1e-6
                xi = [xi,xtmp];
                yi = [yi,ymid];
                Intervals(end+1,:) = [length(xi),Intervals(i,2)];
                Intervals(i,2) = length(xi);
                Flags(end+1) = 1;
            else
                Flags(i) = 0;
            end
        end
        
    end
    [~,index] = sort(xi);
    plot(x,feval(y,x));hold on
    title(num2str(length(xi)))
    plot(sort(xi),yi(index));hold off
    axis off
    pause(0.5)
end

xi = sort(xi);
hold on;
plot(xi,0*xi,'*','color','b')
axis off

