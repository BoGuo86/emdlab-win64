function [] = read_stl_bin(obj,FileDir)

f = fopen(FileDir,'r');

fread(f,80,'*char');

Nt = fread(f,1,'uint32');

n = zeros(3*Nt,1);
p = zeros(9*Nt,1);

for i = 1:Nt
    n(3*i-2:3*i) = fread(f,3,'single');
    p(9*i-8:9*i-6) = fread(f,3,'single');
    p(9*i-5:9*i-3) = fread(f,3,'single');
    p(9*i-2:9*i) = fread(f,3,'single');
    fread(f,1,'uint16');
end

fclose(f);

p = reshape(p,3,[]);
p = p';

obj.facets = 1:3*Nt;
obj.facets = reshape(obj.facets,3,[]);
obj.facets = obj.facets';

[p,~,index] = uniquetol(p,1e-4,'ByRows',true);

obj.facets = index(obj.facets);



end