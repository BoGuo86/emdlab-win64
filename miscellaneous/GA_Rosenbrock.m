
function GA_Rosenbrock

clc;

% initial population
pop = rand(50,2);
c = zeros(4,2);
m = zeros(2,2);

% convergency criteria
errMax = 1e-3;
iterMax = 1000;
err = inf;
iter = 0;

% displaying
close all
p1 = animatedline('color', 'r');
p2 = animatedline('color', 'b');
legend('best fit', 'mean fit');
xlabel('Iterations');
ylabel('Best & Mean');
title('Minimizing Rosenbrock function wit GA');

while iter<iterMax && err>errMax
    
    % generation of new population
    % 1) mating
    c(1:2,:) = mate(pop(1:2,:));
    c(3:4,:) = mate(pop(3:4,:));
    
    % 2) mutation
    for i = 1:2
        m(i,:) = mutation(pop(randi([1,10],1),:));
    end
    
    % new population
    newPop = [pop;c;m];
    
    % eval fitness funcion
    f = fitness(newPop);
    
    % holding strongers
    [f,index] = sort(f);
    pop = newPop(index(1:10),:);
    
    % evaluation of error
    Best = min(f(1:10));
    Mean = mean(f(1:10));
    err = abs(Best - Mean);
    
    % next iteration
    iter = iter + 1;
    
    % displaying
    addpoints(p1, iter, Best);
    addpoints(p2, iter, Mean);
    
end

% best value
disp('best value:');
disp(f(1));

disp('x*:');
disp(pop(1,:));
disp('number of iterations:');
disp(iter);
disp('error:');
disp(err);
box on;

function y = mate(x)
y = x;
y(1,1) = rand*x(1,1) + rand*x(2,1);
y(1,2) = rand*x(1,2) + rand*x(2,2);
y(2,1) = rand*x(1,1) + rand*x(2,1);
y(2,2) = rand*x(1,2) + rand*x(2,2);

function y = mutation(x)
alpha = 2*rand-1;
beta = 2*rand-1;
y = 10*[alpha*x(1), beta*x(2)];

function y = fitness(p)
y = 10*(p(:,2)-p(:,1).^2).^2 + (1-p(:,1)).^2;
