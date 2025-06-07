
function WriteMDataBinary(fdir, mname, material)


f = fopen([fdir,'\',mname,'.dat'],'w');

% material properties
mps = fieldnames(material);

for i = 1:numel(mps)
    mp = mps{i};
    fwrite(f,mp,'char');
    for j = length(mp)+1:80
        fwrite(f,' ','char');
    end
    [m,n] = size(material.(mp));
    fwrite(f,m,'double');
    fwrite(f,n,'double');
    fwrite(f,material.(mp),'double');
end
fclose(f);

end
