
clc
clear

x = 1:10;
y = rand(1,10);


plot(x,y)
hold on

p = pchip(x,y);

x = 1:0.1:10;
plot(x,ppval(p,x))
