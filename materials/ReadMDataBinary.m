function mateial = ReadMDataBinary(fdir,mname)

f = fopen([fdir,'\',mname,'.dat'],'r');

mateial = emdlab_phy_material;
check = true;

while check
    str = fread(f,80,'*char');
    str = strsplit(str',' ');
    if ismember(str{1},{'ThermalConductivity',...
            'ElectricPermitivity',...
            'MagneticPermeability',...
            'MassDensity',...
            'HeatCapacity',...
            'ElectricConductivity',...
            'PoissonRatio',...
            'YoungModulus'})
        
        m = fread(f,1,'double');
        n = fread(f,1,'double');
        value = fread(f,m*n,'double');
        value = reshape(value,m,n);
        mateial.(str{1}) = emdlab_phy_materialProperty(value);
    else
        check = false;
    end
end

fclose(f);

end
