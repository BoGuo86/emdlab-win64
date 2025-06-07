function i_phaseA = SteadyState(linkageFlux,rotorAngle,phaseCurrent,...
    theta_on,theta_off,rpm)

theta_on = theta_on*pi/180;
theta_off = theta_off*pi/180;
w = 2*pi*rpm/60;

di = phaseCurrent(2);
dthetar = rotorAngle(2);

Nthetar = size(linkageFlux,1);
linkageFlux = [linkageFlux;flipud(linkageFlux(2:end,:))];


rotorAngle = [rotorAngle(1:end-1),linspace(pi/6,pi/3,Nthetar)];

Linc = [(linkageFlux(:,2)-linkageFlux(:,1))/di,...
    (linkageFlux(:,3:end)-linkageFlux(:,1:end-2))/2/di,...
    (linkageFlux(:,end)-linkageFlux(:,end-1))/di];

Cw = [(linkageFlux(2,:)-linkageFlux(1,:))/dthetar;...
    (linkageFlux(3:end,:)-linkageFlux(1:end-2,:))/2/dthetar;...
    (linkageFlux(end,:)-linkageFlux(end-1,:))/dthetar];

[rotorAngle,phaseCurrent] = meshgrid(rotorAngle,phaseCurrent);

% parameters
Vdc = 280;
Rdc = 0.6;

h = Hysteron;
h.x_max = 10.5;
h.x_min = 9.5;
h.y_min = 1;
h.y_max = 0;
h.output = 1;

theta = linspace(0,pi/3,2000);
Ntheta = length(theta);
iph = zeros(1,Ntheta);
vph = zeros(1,Ntheta);
lph = zeros(1,Ntheta);

% loop over theta
for i = 2:Ntheta
    
    xLinc = interp2(rotorAngle,phaseCurrent,Linc',theta(i-1),iph(i-1));
    xCw = interp2(rotorAngle,phaseCurrent,Cw',theta(i-1),iph(i-1));
    iph(i) = iph(i-1) + theta(2)*(vph(i-1)-Rdc*iph(i-1)-xCw*w)/(w*xLinc);
    
    % specefying volatge for next iteration
    if theta(i) < theta_on
        vph(i) = 0;
    elseif theta(i) < theta_off

        vph(i) = h.getOutput(iph(i)) * Vdc;
    else
        if iph(i)<0
            iph(i) = 0;
            break
        end
        vph(i) = -Vdc;
    end
    
    lph(i) = interp2(rotorAngle,phaseCurrent,linkageFlux',theta(i),iph(i));
    
end

pout = polyarea(iph,lph)*24*rpm/60;
subplot(221)
plot(theta*180/pi,iph)
xlabel('Rotor Position [Deg]')
set(gca,'xLim',[0 max(theta)*180/pi],'xtick',linspace(0,max(theta)*180/pi,5))
ylabel('Phase Current [A]')
title(['Copper Loss = ',num2str(3*Rdc*sum(iph.^2)*theta(2)*4/pi),' [W]'])
hold on;plot([0,45],[h.x_max,h.x_max],'--');
hold on;plot([0,45],[h.x_min,h.x_min],'--');
subplot(222)
plot(theta*180/pi,vph)
xlabel('Rotor Position [Deg]')
set(gca,'xLim',[0 max(theta)*180/pi],'xtick',linspace(0,max(theta)*180/pi,5))
ylabel('Phase Voltage [V]')
subplot(223)
plot(theta*180/pi,lph)
xlabel('Rotor Position [Deg]')
set(gca,'xLim',[0 max(theta)*180/pi],'xtick',linspace(0,max(theta)*180/pi,5))
ylabel('Linkage Flux [wb]')
subplot(224)
plot(iph,lph)
xlabel('Phase Current [A]')
ylabel('Linkage Flux [wb]')
title(['Rotor Speed = ',num2str(rpm),' [rpm], Output Power = ',num2str(pout),' [W]'])

i_phaseA = [theta;iph];

end
