

clc

f = fopen('D:\My Projects\flat_linear_pm_motor\Maxwell_sims\Normal_Flux.tab', 'r');

if f < 0
    error('can not open file');
end

str = fgetl(f);

x = [];
By = [];

while true
    str = fgetl(f);
    if str<0
        break
    end
    
    str = strsplit(str, '\t');
    x(end+1) = str2num(str{1});
    By(end+1) = str2num(str{2});
    
end
 
fclose(f)


plot(x, By)
    