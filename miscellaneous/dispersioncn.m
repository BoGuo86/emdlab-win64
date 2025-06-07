function result = dispersioncn(q2,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  format long
%   global   f  
i=sqrt(-1)
E2= 2.5
E1= 2.1
i=sqrt(-1)
 pi=4.*atan(1)
 h = 6.6260e-34; %j.S
 hbarr =h./(2.*pi); %j.s
 hh=4.1356e-15% ev.s
 hbar=hh./(2.*pi); %ev.s
 ee=1.6021e-19
 m0=4.*pi.*10.^-7
 c=299792458.*10.^9     % nano m/s
 E0=8.85*1e-12
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 mu=0.5 
  tau=0.5*1e-12
  w=2.*pi.*1e12.*f
%   kT=0.02525;  %%%%%%%%ev
kB=8.6173e-5; %ev/K
T=300; %K
kT=kB*T;

   Gintra=(ee.^2./(hbarr)).*(2.*i.*kT./(pi.*hbar.*(w+i.*tau))).*log(2.*cosh(mu./(2.*kT)));
   Ginter=(ee.^2./(4.*hbarr)).*(1./2+(1./pi).*atan((hbar.*w-2.*mu)./(2.*kT))-(i./(2.*pi)).*...
    log((hbar.*w+2.*mu).^2./(((hbar.*w-2.*mu).^2+(2.*kT).^2))));
sigma=Gintra+Ginter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 150 %nm
V = (w./c).*sqrt(E2- E1)
E0 = 8.36*10.^(6)
Ed = 8.36*10.^(6);
alpha = 10^(-14)

%q1=Sqrt[V^2+q2^2]

P2 = (1/2).*alpha.*E0.^2

ss= sqrt((q2.^4)./(P2.^2) + (4./P2).*(((w./c).*abs(imag(sigma))+ sqrt(V.^2 + q2.^2)).^2 + q2.^2 + P2))

gama= sqrt((ss - (1./P2).*q2.^2)/(2))
delta= sqrt((ss + (1./P2)*q2.^2)/(2))
m= (gama.^2)/(ss)
K=ellipke(m);% calculate complete elliptic integral of the first kind
z0= -sqrt(1./(P2.*ss)).*K
u= sqrt (ss.*P2).*(d + z0)


[s,c,d]=ellipj(u, m);% calculate Jacobian elliptic functions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result=(s.*d)/(c) - ((-w./c).*imag(sigma) + sqrt(V.^2 + q2.^2))/(sqrt(ss.*P2))

end



