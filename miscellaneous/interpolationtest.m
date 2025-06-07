nmname = 'm19_24ga';
f = fopen([cd,'\MaterialsData\bhdata\',nmname,'.txt'],'r');
hb = fscanf(f,'%f');
h = hb(1:2:end);
b = hb(2:2:end);
fclose(f);


plot(h,b,'.')
bb = [b;100];
hh = interp1(b,h,bb,'linear','extrap');
plot(bb,(hh./bb)*4*pi*1e-7)
v = (h./b)*4*pi*1e-7;
vv = interp1(b,v,bb,'linear','extrap');
plot(bb,vv)

