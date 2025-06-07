
function [d,phaseA,phaseB,phaseC] = Winding(Nm,Q)

clc

Nsm = Q/Nm;
Ncph = Q/3;

for K0 = 1:Q
    if rem(3*Nm*K0/2/Q,3) == 1
        break
    end
    if K0 == Q
        error('No phase offset can be found.');
    end
end

S = max(floor(Nsm),1);
slot = (0:Q-1)';
d = ones(Q,1);
In = slot;
Out = circshift(slot,-S);
Angle = (slot*Nm/Q)*180;
Angle = rem(Angle+180,360)-180;
index = Angle>90;
Angle(index) = Angle(index)-180;
d(index) = -1;
index = Angle<-90;
Angle(index) = Angle(index)+180;
d(index) = -1;

fprintf('Phase Shift = %d\n',K0)
fprintf('Coil Span = %d\n',S)

% [~,index] = sort(abs(Angle));
clc
disp([In,Angle])
phaseA = input('klj');
% phaseA = index(1:Ncph);

d = d(phaseA);
phaseA = [In(phaseA),Out(phaseA)];

phaseB = rem(phaseA+K0,Q);
phaseC = rem(phaseB+K0,Q);


disp(table(phaseA,phaseB,phaseC,...
    'VariableNames',{'phaseA','phaseB','phaseC'}));

hold all
for i = 1:size(phaseA,1)
    x = cos(phaseA(i,1)*2*pi/Q);
    y = sin(phaseA(i,1)*2*pi/Q);
    plot(x,y);
    if d(i)>0
        text(x,y,'A+','color','b');
    else
        text(x,y,'A-','color','b');
    end
    
    x = 1.2*cos(phaseA(i,2)*2*pi/Q);
    y = 1.2*sin(phaseA(i,2)*2*pi/Q);
    plot(x,y);
    if d(i)>0
        text(x,y,'A-','color','b');
    else
        text(x,y,'A+','color','b');
    end
    
    x = cos(phaseB(i,1)*2*pi/Q);
    y = sin(phaseB(i,1)*2*pi/Q);
    plot(x,y);
    if d(i)>0
        text(x,y,'B+','color','r');
    else
        text(x,y,'B-','color','r');
    end
    
    x = 1.2*cos(phaseB(i,2)*2*pi/Q);
    y = 1.2*sin(phaseB(i,2)*2*pi/Q);
    plot(x,y);
    if d(i)>0
        text(x,y,'B-','color','r');
    else
        text(x,y,'B+','color','r');
    end
    
    x = cos(phaseC(i,1)*2*pi/Q);
    y = sin(phaseC(i,1)*2*pi/Q);
    plot(x,y);
    if d(i)>0
        text(x,y,'C+','color','m');
    else
        text(x,y,'C-','color','m');
    end
    
    x = 1.2*cos(phaseC(i,2)*2*pi/Q);
    y = 1.2*sin(phaseC(i,2)*2*pi/Q);
    plot(x,y);
    if d(i)>0
        text(x,y,'C-','color','m');
    else
        text(x,y,'C+','color','m');
    end
end
axis off equal

end

