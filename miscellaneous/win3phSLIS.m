% single layer integral slot three phase circuit
function [phaseA,phaseB,phaseC] = win3phSLIS(p,Q)

if rem(p,2)
    error('number of poles must be an even number.');
end

% number of phases
m = 3;
% number of slot per pole per phase
q = Q/p/m;
% checking for integral slot
if (q-floor(q))
    error('q is not an integer, in integral slot windings q must be an integer.');
end

% span
S = m*q;

slot = (0:Q-1)';

slot = reshape(slot,S,[]);

phaseA_slot = slot(1:q,:);
phaseA_in = phaseA_slot(:,1:2:end);
phaseA_out = phaseA_slot(:,2:2:end);
phaseA = [phaseA_in(:),phaseA_out(:)];

K0 = 2*Q/3/p;

phaseB = rem(phaseA+K0,Q);
phaseC = rem(phaseB+K0,Q);

clc
disp(table(phaseA,phaseB,phaseC,...
    'VariableNames',{'phaseA','phaseB','phaseC'}));

hold all
for i = 1:size(phaseA,1)
    x = cos(phaseA(i,1)*2*pi/Q);
    y = sin(phaseA(i,1)*2*pi/Q);
    plot(x,y);
    text(x,y,'A+','color','b');
    x = cos(phaseA(i,2)*2*pi/Q);
    y = sin(phaseA(i,2)*2*pi/Q);
    plot(x,y);
    text(x,y,'A-','color','b');
    x = cos(phaseB(i,1)*2*pi/Q);
    y = sin(phaseB(i,1)*2*pi/Q);
    plot(x,y);
    text(x,y,'B+','color','r');
    x = cos(phaseB(i,2)*2*pi/Q);
    y = sin(phaseB(i,2)*2*pi/Q);
    plot(x,y);
    text(x,y,'B-','color','r');
    x = cos(phaseC(i,1)*2*pi/Q);
    y = sin(phaseC(i,1)*2*pi/Q);
    plot(x,y);
    text(x,y,'C+','color','m');
    x = cos(phaseC(i,2)*2*pi/Q);
    y = sin(phaseC(i,2)*2*pi/Q);
    plot(x,y);
    text(x,y,'C-','color','m');
end
axis off equal

title(['[3ph Integral Slot Winding]','[Q = ',num2str(Q),']',...
    '[p = ',num2str(p),']'])

end



