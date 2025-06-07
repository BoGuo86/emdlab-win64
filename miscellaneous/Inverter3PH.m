clc
clear

t = 0:1e-4:0.1;
Nt = length(t);
w = 2*pi*50;

va = zeros(1,Nt);
vb = zeros(1,Nt);
vc = zeros(1,Nt);


pa = pulse_func();
pb = pulse_func();
pb.t0 = 0.02/3;
pc = pulse_func();
pc.t0 = 2*0.02/3;

Vdc = 400;
VdcP2 = Vdc/2;

for i = 1:Nt
    s1 = pa.getValue(t(i));
    s2 = pb.getValue(t(i));
    s3 = pc.getValue(t(i));
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
    else
        error('Short Circuit');
    end
    va(i) = (2*va0-vb0-vc0)/3;
    vb(i) = (-va0+2*vb0-vc0)/3;
    vc(i) = (-va0-vb0+2*vc0)/3;
end

Resistance = 1;
Inductance = 1e-3;

subplot(311)
plot(t,va,t,rl_vfed(va,Resistance,Inductance,t(2)))
subplot(312)
plot(t,vb,t,rl_vfed(vb,Resistance,Inductance,t(2)))
subplot(313)
plot(t,vc,t,rl_vfed(vc,Resistance,Inductance,t(2)))