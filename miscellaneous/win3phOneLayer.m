
function [In,Out,phaseA,phaseB,phaseC] = win3phOneLayer(p,Q)

if rem(p,2)
    error('number of poles must be an even number.');
end


% number of phases
m = 3;
% number of slot per pole per phase
q = Q/p/m;
% checking for integral slot
if ~isinteger(uint8(q))
    error('q is not an integer, in integral slot windings q must be an integer.');
end

% span
S = m*q;

Coil = (1:Q/2)';

In = zeros(length(Coil),1);

for i = 1:p/2
    In(1+(i-1)*S:i*S) = (1:S)+2*(i-1)*S;
end

Out = In + S;

Angle = ((Coil-1)*p/Q)*360;
Angle = rem(Angle,360);

Angle = rem(Angle,180);
[~,index] = sort(Angle);

% fprintf('Phase Shift = %d\n',K0)
fprintf('Coil Span = %d\n',S)
disp(table(Coil,In,Out,Angle,...
    'VariableNames',{'Coil','In','Out','Angle'}));

Ncph = Q/2/m;
phaseA = index(1:Ncph);
phaseB = phaseA + Ncph ;
phaseC = phaseB + Ncph ;

index = true;
while index
    index = phaseB>Q/2;
    phaseB(index) = phaseB(index)-Q/2;
end


index = true;
while index
    index = phaseC>Q/2;
    phaseC(index) = phaseC(index)-Q/2;
end

disp(table(phaseA,phaseB,phaseC,...
    'VariableNames',{'phaseA','phaseB','phaseC'}));






end