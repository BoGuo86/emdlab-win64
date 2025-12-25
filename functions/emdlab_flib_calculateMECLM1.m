function [totalLoss,pevec] = emdlab_flib_calculateMECLM1(obj, plotFlag)

if nargin<2, plotFlag = false; end

pevec = 0;
totalLoss = 0;
for mzName = string(fieldnames(obj.m.mzs)')

    mzptr = obj.m.mzs.(mzName);
    if mzptr.props.isEddyZone

        % get simulation time
        dt = obj.simTime(2);
        gv_EPeriod = obj.simTime(end);

        Az = mzptr.props.Az;
        dAzdt = zeros(size(Az));

%         for i = 1:size(Az,1)
%             pp  = pchip(obj.simTime, Az(i,:));   % shape-preserving cubic
%             ppd = fnder(pp, 1);       % derivative
%             dAzdt(i,:) = ppval(ppd, obj.simTime);
%         end
%         dAzdt(:,1) = 0;

        dAzdt = [0*Az(:,1),...
            (Az(:,3:end)-Az(:,1:end-2))/2,...
            Az(:,end)-Az(:,end-1)]/dt;

        e_mvolume = (mzptr.getAreaOfElements*obj.getDepth*1e-6);
        peddy = zeros(mzptr.Ne,length(obj.simTime));
        for i = 1:length(obj.simTime)
            dAzdt_i = dAzdt(:,i);
            dAzdt_i = dAzdt_i(mzptr.cl);
            peddy(:,i) = (1/6)* (sum(dAzdt_i.^2,2) + ...
                dAzdt_i(:,1).*dAzdt_i(:,2)+dAzdt_i(:,2).*dAzdt_i(:,3)+dAzdt_i(:,3).*dAzdt_i(:,1));
        end

        peddy = peddy * obj.m.mts.(mzptr.material).ElectricConductivity.value * dt .* e_mvolume/ gv_EPeriod;
        pevec = pevec + sum(peddy);
        peddy = sum(peddy,2);

        mzptr.props.peddy = sum(peddy);
        totalLoss = totalLoss + mzptr.props.peddy;

        disp(mzName + " --> Eddy Current Loss = " + sprintf('%.2f', mzptr.props.peddy ));

        if plotFlag
            obj.plotCoreLossDensity(mzName,(peddy./e_mvolume)');
            title("Eddy Current Loss Density of <" + mzName +">")
        end

    end

end

pevec = pevec*gv_EPeriod/dt;
disp("*** Total Loss = " + sprintf('%.2f', totalLoss) + " ***");

end

