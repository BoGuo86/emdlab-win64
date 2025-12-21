function emdlab_flib_calculateIronLossesModel2(obj, Kh, alpha, beta, Ke, plotFlag)

if nargin<6, plotFlag = true; end

totalCoreLoss = 0;
for mzName = string(fieldnames(obj.m.mzs)')

    mzptr = obj.m.mzs.(mzName);
    if mzptr.props.isCoreLossActivated

        % calculate hysteresis & excess losses
        % get required fields data
        bx = mzptr.props.Bxg';
        by = mzptr.props.Byg';

        % get simulation time
        dt = obj.simTime(2);
        gv_EPeriod = obj.simTime(end);

        % get mass density of iron material
        mass_density = obj.m.mts.(mzptr.material).MassDensity.value;

        % calculation of KGSE
        temp = linspace(0,2*pi,1000);
        KGSE = Kh/((2*pi)^(alpha-1))/(sum((abs(cos(temp))).^alpha .* ...
            (abs(sin(temp))).^(beta-alpha))*temp(2));

        bmag = sqrt(bx.^2+by.^2);
        dbmagdt = [bmag(:,2)-bmag(:,1),...
            (bmag(:,3:end)-bmag(:,1:end-2))/2,...
            bmag(:,end)-bmag(:,end-1)]/dt;

        dbxdt = [bx(:,2)-bx(:,1),...
            (bx(:,3:end)-bx(:,1:end-2))/2,...
            bx(:,end)-bx(:,end-1)]/dt;

        dbydt = [by(:,2)-by(:,1),...
            (by(:,3:end)-by(:,1:end-2))/2,...
            by(:,end)-by(:,end-1)]/dt;

        phystx = (abs(dbmagdt).^alpha) .* (bmag.^(beta-alpha));
        physty = (abs(dbmagdt).^alpha) .* (bmag.^(beta-alpha));

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

