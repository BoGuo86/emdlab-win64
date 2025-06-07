clc
clear

x = 0:255;

sa = sin(x/255 * 2*pi);
sb = sin(x/255 * 2*pi - 2*pi/3);
sc = sin(x/255 * 2*pi - 4*pi/3);

sa = (sa+1)/2; sa = floor(sa*255);
sb = (sb+1)/2; sb = floor(sb*255);
sc = (sc+1)/2; sc = floor(sc*255);

plot(x,sa,x,sb,x,sc)

clc
index = 0;
for i = 1:16
    for j = 1:16
        index = index + 1;
        fprintf('%d,',sb(index));
    end
    fprintf('\n');
end
legend({'a','b','c'})