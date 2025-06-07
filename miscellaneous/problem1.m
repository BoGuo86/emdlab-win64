clc
clear
t = [1;2;3;4;5;6];
y = [0.5;2.2;5.3;6.1;7;9.4];
hold all
plot(t,y,'*');
c = fit(t,y,'poly1');
% plot(c);
c = fit(t,y,'poly2');
% plot(c);
c = fit(t,y,'poly3');
% plot(c);
