
t = 0:1e-5:0.1;

va = t;
vb = t;
vc = t;

ca = 0*t;
cb = 0*t;
cc = 0*t;

ea = xy_func(rotorAngle*p.wm,BEMFs(1,:));
eb = xy_func(rotorAngle*p.wm,BEMFs(2,:));
ec = xy_func(rotorAngle*p.wm,BEMFs(3,:));

s = tri_func(-1,1,1e4);
fa = cos_func(1,97.5,0);
fb = cos_func(1,97.5,-2*pi/3);
fc = cos_func(1,97.5,-4*pi/3);

Vdc = 400;




R = 8;
L = 0.090464494127805;

tmp1 = L/t(2);
tmp2 = L/t(2)+R;


for i = 2:length(t)
    
    b1 = fa.getValue(t(i))>s.getValue(t(i));
    b2 = fb.getValue(t(i))>s.getValue(t(i));
    b3 = fc.getValue(t(i))>s.getValue(t(i));
    
    if b1&&~b2
        va(i) = Vdc/2;vb(i) = -Vdc/2;vc(i) = 0;
    elseif b1&&~b3
        va(i) = Vdc/2;vb(i) = 0;vc(i) = -Vdc/2;
    elseif b2&&~b1
        va(i) = -Vdc/2;vb(i) = Vdc/2;vc(i) = 0;
    elseif b2&&~b3
        va(i) = 0;vb(i) = Vdc/2;vc(i) = -Vdc/2;
    elseif b3&&~b1
        va(i) = -Vdc/2;vb(i) = 0;vc(i) = Vdc/2;
    else
        va(i) = 0;vb(i) = -Vdc/2;vc(i) = Vdc/2;
    end
    
    ca(i) = (va(i)-ea.getValue(t(i))+ca(i-1)*tmp1)/tmp2;
    cb(i) = (vb(i)-eb.getValue(t(i))+cb(i-1)*tmp1)/tmp2;
    cc(i) = (vc(i)-ec.getValue(t(i))+cc(i-1)*tmp1)/tmp2;
    
end


plot(t,ca,t,cb,t,cc)