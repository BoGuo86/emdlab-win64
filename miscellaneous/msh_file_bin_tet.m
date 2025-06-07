clear
clc

f = fopen('C:\Users\Ali Jamalifard\Desktop\ali.msh','r');

if f<0
    error('can not open file');
end

e = cell(1,200);
for i = 1:numel(e)
e{i} = zeros(4,[]);
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
                    fread(f,Nt,'int32');
                    fread(f,3,'int32');
                end
            elseif eIndex == 4
                Ne = fread(f,1,'int32');
                Nt = fread(f,1,'int32');
                for j = 1:Ne
                    fread(f,1,'int32');
                    fread(f,Nt-1,'int32');
                    zIndex = fread(f,1,'int32');
                    e{zIndex}(:,end+1) = fread(f,4,'int32');
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
p = p';
tic
m = TTMDBC;
for i = 1:zIndex
    xtmp = p(:,1);
    ytmp = p(:,2);
    ztmp = p(:,3);
    
    xtmp = xtmp(e{i});
    ytmp = ytmp(e{i});
     ztmp = ztmp(e{i});
     
    [ptmp, ~, index] = uniquetol([xtmp(:),ytmp(:),ztmp(:)], 1e-6, 'ByRows', true);
    
    etmp = 1:4*size(e{i},2);
    etmp = etmp(index);
    etmp = reshape(etmp,4,[])';
    
    m.addmz(['z',num2str(i)],TTMZPC(etmp,ptmp));
end
m.showmzs
toc
