clear all
clc

f = fopen('C:\Users\AliJamalifard\Desktop\gmesh\INduction.msh','r');

e = cell(1,10);
for i = 1:numel(e)
e{i} = zeros(3,[]);
end

while true
    str =  fgetl(f);
    if strcmp(str,'$Nodes')
        Nn = fgetl(f);
        Nn = str2double(Nn);
        p = zeros(3,Nn);
        for i = 1:Nn
            fread(f,1,'int32');
            p(:,i) = fread(f,3,'double');
        end
    end
    if strcmp(str,'$Elements')
        Nn = fgetl(f);
        Nn = str2double(Nn);
        for i = 1:Nn
            eIndex = fread(f,1,'int32');
            if eIndex == 15
                Ne = fread(f,1,'int32');
                Nt = fread(f,1,'int32');
                for j = 1:Ne
                    fread(f,1,'int32');
                    fread(f,Nt,'int32');
                    fread(f,1,'int32');
                end
            elseif eIndex == 1
                Ne = fread(f,1,'int32');
                Nt = fread(f,1,'int32');
                for j = 1:Ne
                    fread(f,1,'int32');
                    fread(f,Nt,'int32');
                    fread(f,2,'int32');
                end
            elseif eIndex == 2
                Ne = fread(f,1,'int32');
                Nt = fread(f,1,'int32');
                for j = 1:Ne
                    fread(f,1,'int32');
                    fread(f,Nt-1,'int32');
                    zIndex = fread(f,1,'int32');
                    e{zIndex}(:,end+1) = fread(f,3,'int32');
                end
            end
        end
    end
    if strcmp(str,'$EndElements')
        break;
    end
end
fclose(f);
e = e(1:zIndex);

p = p(1:2,:)';

tic
m = TMDBC;
for i = 1:zIndex
    xtmp = p(:,1);
    ytmp = p(:,2);
    xtmp = xtmp(e{i});
    ytmp = ytmp(e{i});
    [ptmp, ~, index] = uniquetol([xtmp(:),ytmp(:)], 1e-6, 'ByRows', true);
    
    etmp = 1:3*size(e{i},2);
    etmp = etmp(index);
    etmp = reshape(etmp,3,[])';
    
    m.addmz(TMZPC(etmp,ptmp));
end
m.showmzs
toc
% m = TMDBC;
% m.addmz('s',TMZPC(el(1).e,p));
% m.addmz('s1',TMZPC(el(2).e,p));
