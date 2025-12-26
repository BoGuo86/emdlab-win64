function totalCoreLoss = emdlab_flib_calculateIronLossesModel1(obj, Kh, alpha, beta, Ke, plotFlag)

if nargin<6, plotFlag = false; end

totalCoreLoss = 0;
for mzName = string(fieldnames(obj.m.mzs)')

    mzptr = obj.m.mzs.(mzName);
    if mzptr.props.isCoreLossActivated

        % calculate hysteresis & excess losses
        % get required fields data
        bx = mzptr.props.Bxg';
        by = mzptr.props.Byg';

        

        if mzptr.props.isMoving

            
            xsh = [0;obj.movingRegions.(mzptr.props.movingRegionName).motionHistory(:,1)];
            ysh = [0;obj.movingRegions.(mzptr.props.movingRegionName).motionHistory(:,2)];
            xc = [0;obj.movingRegions.(mzptr.props.movingRegionName).motionHistory(:,3)];
            yc = [0;obj.movingRegions.(mzptr.props.movingRegionName).motionHistory(:,4)];
            theta = [0;obj.movingRegions.(mzptr.props.movingRegionName).motionHistory(:,5)];

            xaxis = [1;0];
            bx_new = bx;
            by_new = by;
            for i = 2:size(bx,2)
                xaxis = xaxis + [xsh(i);ysh(i)];
                xaxis = [xc(i);yc(i)] + [cos(theta(i)),-sin(theta(i));sin(theta(i)),cos(theta(i))] * (xaxis-[xc(i);yc(i)]);
                yaxis = [-xaxis(2);xaxis(1)];
                bx_new(:,i) = bx(:,i)*xaxis(1) + by(:,i)*xaxis(2);
                by_new(:,i) = bx(:,i)*yaxis(1) + by(:,i)*yaxis(2);
            end

            theta = cumsum(theta);
bx_new = bx.*cos(theta') + by.*sin(theta');
by_new = -bx.*sin(theta') + by.*cos(theta');

            bx = bx_new;
            by = by_new;

        end

        bx = bx-mean(bx,2);
            by = by-mean(by,2);

        % get simulation time
        dt = obj.simTime(2);
        gv_EPeriod = obj.simTime(end);

        % get mass density of iron material
        mass_density = obj.m.mts.(mzptr.material).MassDensity.value;

        % calculation of KGSE
        temp = linspace(0,2*pi,1000);
        KGSE = Kh/((2*pi)^(alpha-1))/(sum((abs(cos(temp))).^alpha .* ...
            (abs(sin(temp))).^(beta-alpha))*temp(2));

        dbxdt = [bx(:,2)-bx(:,1),...
            (bx(:,3:end)-bx(:,1:end-2))/2,...
            bx(:,end)-bx(:,end-1)]/dt;

        dbydt = [by(:,2)-by(:,1),...
            (by(:,3:end)-by(:,1:end-2))/2,...
            by(:,end)-by(:,end-1)]/dt;

        phystx = (abs(dbxdt).^alpha) .* (abs(bx).^(beta-alpha));
        physty = (abs(dbydt).^alpha) .* (abs(by).^(beta-alpha));

        phystx(:,1) = phystx(:,1)/2;
        phystx(:,end) = phystx(:,end)/2;

        physty(:,1) = physty(:,1)/2;
        physty(:,end) = physty(:,end)/2;

        phystx = sum(phystx ,2);
        physty = sum(physty ,2);
        physt = phystx+physty;

        e_mass = mass_density * (mzptr.getAreaOfElements*obj.getDepth*1e-6);
        physt = physt * KGSE * dt .* e_mass/ gv_EPeriod;

        mzptr.props.physt = sum(physt);
        totalCoreLoss = totalCoreLoss + mzptr.props.physt;

        disp(mzName + " --> Hysteresis Loss = " + sprintf('%.2f', mzptr.props.physt));

        if plotFlag
            obj.plotCoreLossDensity(mzName,(physt*mass_density./e_mass)');
            title("Hysteresis & Excess Loss Density of <" + mzName +">")
        end

        % calculate eddy current losses
        peddy = (dbxdt.^2+dbydt.^2);

        peddy(:,1) = peddy(:,1) / 2;
        peddy(:,end) = peddy(:,end) / 2;

        peddy = sum(peddy,2);
        e_mass = mass_density * (mzptr.getAreaOfElements*obj.getDepth*1e-6);

        peddy =  (Ke/2/pi^2)*peddy * dt .* e_mass/ gv_EPeriod;

        mzptr.props.peddy = sum(peddy);
        totalCoreLoss = totalCoreLoss + mzptr.props.peddy;

        disp(mzName + " --> Eddy Current Loss = " + sprintf('%.2f', mzptr.props.peddy));

        if plotFlag
            obj.plotCoreLossDensity(mzName,(peddy*mass_density./e_mass)');
            title("Eddy Current Loss Density of <" + mzName +">")
        end

    end

end

disp("*** Total Core Loss = " + sprintf('%.2f', totalCoreLoss) + " ***");

end

