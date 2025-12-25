function totalLoss = emdlab_flib_calculateMECLM2(obj, plotFlag)

if nargin<2, plotFlag = false; end

totalLoss = 0;

t  = obj.simTime;
Nt = length(t);
T  = t(end);

% 2-point Gauss rule (reference interval [-1,1])
gauss_x = [-1/sqrt(3), 1/sqrt(3)];
gauss_w = [1, 1];

for mzName = string(fieldnames(obj.m.mzs)')

    mzptr = obj.m.mzs.(mzName);

    if mzptr.props.isEddyZone

        Az = mzptr.props.Az;    % nodal A_z, size: Nnodes x Nt
        sigma = obj.m.mts.(mzptr.material).ElectricConductivity.value;

        % build pchip derivative for each node (once)
        ppd = cell(size(Az,1),1);
        for i = 1:size(Az,1)
            pp  = spline(t, Az(i,:));
            ppd{i} = fnder(pp,1);   % quadratic in time
        end

        elemArea = mzptr.getAreaOfElements;     % [mm^2]
        depth    = obj.getDepth;                % [mm]
        e_mvolume = elemArea .* depth * 1e-6;   % [m^3]

        elemEnergy = zeros(mzptr.Ne,1);

        % --- time integration ---
        for k = 1:Nt-1

            tk  = t(k);
            tk1 = t(k+1);
            dt  = tk1 - tk;

            % map Gauss points to [tk, tk1]
            tg = 0.5*(tk1-tk)*gauss_x + 0.5*(tk1+tk);

            for g = 1:2

                dAzdt_g = zeros(size(Az,1),1);
                for i = 1:size(Az,1)
                    dAzdt_g(i) = ppval(ppd{i}, tg(g));
                end

                % element-wise nodal derivatives
                dAzdt_e = dAzdt_g(mzptr.cl);   % Ne x 3

                % exact spatial integration for linear triangle
                p_e = (1/6) * ( ...
                    sum(dAzdt_e.^2,2) + ...
                    dAzdt_e(:,1).*dAzdt_e(:,2) + ...
                    dAzdt_e(:,2).*dAzdt_e(:,3) + ...
                    dAzdt_e(:,3).*dAzdt_e(:,1) );

                % accumulate energy
                elemEnergy = elemEnergy + ...
                    gauss_w(g) * p_e * (dt/2);
            end
        end

        % multiply physical factors
        elemEnergy = sigma .* elemEnergy .* e_mvolume;

        % average power
        peddy_elem = elemEnergy / T;
        mzptr.props.peddy = sum(peddy_elem);
        totalLoss = totalLoss + mzptr.props.peddy;

        disp(mzName + " --> Eddy Current Loss = " + ...
             sprintf('%.3f', mzptr.props.peddy));

        if plotFlag
            obj.plotCoreLossDensity(mzName,(peddy_elem./e_mvolume)');
            title("Eddy Current Loss Density of <" + mzName + ">");
        end

    end
end

disp("*** Total Loss = " + sprintf('%.3f', totalLoss) + " ***");

end
