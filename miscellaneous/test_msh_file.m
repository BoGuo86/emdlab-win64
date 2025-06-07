
f = fopen('C:\Users\Ali Jamalifard\Desktop\geom.msh','r');
while true
    str =  fgetl(f);
    if strcmp(str,'$Nodes')
        Nn = fgetl(f);
        Nn = str2num(Nn);
        p = zeros(Nn,3);
        for i = 1:Nn
            fread(f,1,'int32');
            p(i,:) = fread(f,3,'double');
        end
    elseif strcmp(str,'$Elements')
        Nn = fgetl(f);
        Nn = str2num(Nn);
        e = zeros(3,[]);
        a = zeros(2,[]);
        for i = 1:Nn
            
            
            if fread(f,1,'int32') == 15
                Ne = fread(f,1,'int32');
                Nt = fread(f,1,'int32');
                for j = 1:Ne
                    fread(f,1,'int32');
                    fread(f,Nt,'int32');
                    fread(f,1,'int32');
                end
            elseif fread(f,1,'int32') == 1
                Ne = fread(f,1,'int32');
                Nt = fread(f,1,'int32'); 
                for j = 1:Ne
                    fread(f,1,'int32');
                    fread(f,Nt,'int32');
                    a = [a,fread(f,2,'int32')];
                end
            elseif fread(f,1,'int32') == 2
                disp('sdg')
                Ne = fread(f,1,'int32');
                Nt = fread(f,1,'int32');
                for j = 1:Ne
                    fread(f,1,'int32');
                    fread(f,Nt,'int32');
                    e = [e, fread(f,3,'int32')];
                end      
            end
        end
    elseif strcmp(str,'$EndElements')        
        break;
    end 
end
fclose(f);
