clear all
clc

f = fopen('C:\Users\AliJamalifard\Desktop\gmesh\tmp.msh','r');

el(2) = struct('e', zeros([],3));

while true
    str =  fgetl(f);
    if strcmp(str,'$Nodes')
        Nn = fgetl(f);
        Nn = str2double(Nn);
        p = zeros(Nn,2);
        for i = 1:Nn
            str =  strsplit(fgetl(f), ' ');
            p(i,:) = [str2double(str{2}),str2double(str{3})];
        end
    elseif strcmp(str,'$Elements')
        Nn = fgetl(f);
        Nn = str2double(Nn);
%         e = zeros([], 3);
        
        for i = 1:Nn
            str =  strsplit(fgetl(f), ' ');
            index = str2double(str{end-3});
            if strcmp(str{2},'2')
                el(index).e = [el(index).e;str2double(str{end-2}),str2double(str{end-1}),str2double(str{end})];
            end
        end
    elseif strcmp(str,'$EndElements')
        break;
    end
end
fclose(f);
m = TMDBC;
m.addmz('s',TMZPC(el(1).e,p));
m.addmz('s1',TMZPC(el(2).e,p));
