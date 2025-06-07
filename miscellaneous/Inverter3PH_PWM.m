clc
clear

t = 0:1e-5:0.1;
Nt = length(t);
w = 2*pi*50;

va = zeros(1,Nt);
vb = zeros(1,Nt);
vc = zeros(1,Nt);


pa = cos_func();
pa.A = 1;
pa.f = 25;
pb = pa;
pb.phi0 = -2*pi/3;
pc = pa;
pc.phi0 = -4*pi/3;

f = tri_func();
f.f = 1e3;

Vdc = 400;
VdcP2 = Vdc/2;

for i = 1:Nt
    tmp = f.getValue(t(i));
    s1 = pa.getValue(t(i))>tmp;
    s2 = pb.getValue(t(i))>tmp;
    s3 = pc.getValue(t(i))>tmp;
    if s1&&~s2&&~s3
        va0 = VdcP2; vb0 = -VdcP2; vc0 = -VdcP2;
    elseif s1&&s2&&~s3
        va0 = VdcP2; vb0 = VdcP2; vc0 = -VdcP2;
    elseif s1&&~s2&&s3
        va0 = VdcP2; vb0 = -VdcP2; vc0 = VdcP2;
    elseif ~s1&&s2&&~s3
        va0 = -VdcP2; vb0 = VdcP2; vc0 = -VdcP2;
    elseif ~s1&&s2&&s3
        va0 = -VdcP2; vb0 = VdcP2; vc0 = VdcP2;
    elseif ~s1&&~s2&&s3
        va0 = -VdcP2; vb0 = -VdcP2; vc0 = VdcP2;
    end
    va(i) = (2*va0-vb0-vc0)/3;
    vb(i) = (-va0+2*vb0-vc0)/3;
    vc(i) = (-va0-vb0+2*vc0)/3;
end

Resistance = 8;
Inductance = 80e-3;

subplot(311)
ia = rl_vfed(va,Resistance,Inductance,t(2));
Vm = max(va);Im = max(ia);
plot(t,va/Vm,t,ia/Im);set(gca,'ylim',[-1.5,1.5])
title(['vA_{max} = ',num2str(Vm),' ,iA_{max} = ',num2str(Im)])

subplot(312)
ib = rl_vfed(vb,Resistance,Inductance,t(2));
Vm = max(vb);Im = max(ib);
plot(t,vb/Vm,t,ib/Im);set(gca,'ylim',[-1.5,1.5])
title(['vB_{max} = ',num2str(Vm),' ,iB_{max} = ',num2str(Im)])

subplot(313)
ic = rl_vfed(vc,Resistance,Inductance,t(2));
Vm = max(vc);Im = max(ic);
plot(t,vc/Vm,t,ic/Im);set(gca,'ylim',[-1.5,1.5])
title(['vC_{max} = ',num2str(Vm),' ,iC_{max} = ',num2str(Im)])