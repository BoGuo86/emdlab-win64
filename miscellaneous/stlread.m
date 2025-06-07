
f = fopen('C:\Users\AliJamalifard\Desktop\Part3.STL','r');

n = zeros([],3);
p = zeros([],3);

while true
    str = fgets(f);
    if isnumeric(str)
        break
    end
    str = strsplit(str,' ');
    if strcmp(str{1},'solid')
        while true
            str = strsplit(fgets(f),' ');
            if strcmp(str{1},'endsolid')
                break
            end
            n = [n;str2double(str{4}),str2double(str{5}),str2double(str{6})];
            fgets(f);
            str = strsplit(fgets(f),' ');
            p = [p;str2double(str{3}),str2double(str{4}),str2double(str{5})];
            str = strsplit(fgets(f),' ');
            p = [p;str2double(str{3}),str2double(str{4}),str2double(str{5})];
            str = strsplit(fgets(f),' ');
            p = [p;str2double(str{3}),str2double(str{4}),str2double(str{5})];
            fgets(f);
            fgets(f);      
        end
    end
end

fclose(f);

t = 1:size(p,1);
t = reshape(t,3,[]);
t = t';

[p,~,index] = uniquetol(p,1e-4,'ByRows',true);

trisurf(index(t),p(:,1),p(:,2),p(:,3),'FaceColor','b','edgeColor','w');

set(gcf,'Color','k')
 set(gcf,'Renderer','opengl');
axis equal off
rotate3d on

DT = delaunayTriangulation(p);
faceColor  = [0.6875 0.8750 0.8984];
tetramesh(DT,'FaceColor',faceColor,'FaceAlpha',0.3);
axis off equal