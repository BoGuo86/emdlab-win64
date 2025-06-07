

x = 0:0.01:pi;
y = sin(x);

tic
dy = int_trap(y,x);
toc

