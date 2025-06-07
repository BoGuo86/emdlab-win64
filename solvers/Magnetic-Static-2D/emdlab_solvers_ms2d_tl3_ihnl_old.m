% developer: https://ComProgExpert.com, Ali Jamali-Fard
% two dimensional nonlinear magnetic-static solver
% nonlinear solver: Newton-Raphson
% first order triangular mesh
% triangular lagrangian elements: 3 points per element
% isotropic
% homogenous
% nonlinear

classdef emdlab_solvers_ms2d_tl3_ihnl_old < handle & emdlab_solvers_ms2d_tl3
    
    methods
        
        % Initialization
        function obj = emdlab_solvers_ms2d_tl3_ihnl_old(m)
            
            % mesh pointer
            obj.m = m;
            
            % default settings for solver
            obj.solverSettings.relativeError = 1e-4;
            obj.solverSettings.maxIteration = 20;
            
            % set default properties of mzs
            mzNames = fieldnames(obj.m.mzs);
            
            for i = 1:numel(mzNames)
                obj.setdp(mzNames{i});
            end
            
        end
        
        % Solver
        function setSolverMaxIteration(obj, maxIteration)
            
            if maxIteration < 0 || rem(maxIteration, 1)
                error('maxIteration must be a positive integer.');
            end
            
            obj.solverSettings.maxIteration = maxIteration;
            
        end
        
        function setSolverRelativeError(obj, relativeError)
            
            if relativeError < 0
                error('relativeError must be a real positive number.');
            end
            
            obj.solverSettings.relativeError = relativeError;
            
        end
        
        function setMonitor(obj, value)
            
            obj.monitorResiduals = value;
            
        end
        
        function assignEdata(obj, InitNur)
            
            % check states
            if obj.isElementDataAssigned, return; end
            
            % preparing mesh data
            obj.m.evalKexy1_TL3;
            tic, disp('-------------------------------------------------------');
            
            % assigning material and force data to each triangle            
            % allocation of memory
            obj.edata.MagneticReluctivity = zeros(1, obj.m.Ne);
            obj.edata.ElectricConductivity = zeros(1, obj.m.Ne);
            obj.edata.InternalCurrentDensity = zeros(1, obj.m.Ne);
            obj.edata.MagnetizationX = zeros(1, obj.m.Ne);
            obj.edata.MagnetizationY = zeros(1, obj.m.Ne);
            
            % getting mesh zones
            mzsName = fieldnames(obj.m.mzs);
            
            % loop over mesh zones
            for i = 1:obj.m.Nmzs
                
                mzptr = obj.m.mzs.(mzsName{i});
                
                if ~obj.m.mts.(mzptr.material).MagneticPermeability.isIsotropic
                    
                    throw(MException('', 'Some materials are non-isotropic, please select the correct solver.'));
                    
                elseif obj.m.mts.(mzptr.material).MagneticPermeability.isLinear
                    
                    % assigning Magnetic Permeability
                    obj.edata.MagneticReluctivity(obj.m.ezi(:, mzptr.zi)) = 1/obj.m.mts.(mzptr.material).MagneticPermeability.value;
                    
                else
                    
                    if nargin == 2
                        obj.edata.MagneticReluctivity(obj.m.ezi(:, mzptr.zi)) = InitNur * obj.pcts.nu0;
                    else
                        obj.edata.MagneticReluctivity(obj.m.ezi(:, mzptr.zi)) = 0.001 * obj.pcts.nu0;
                    end
                    
                end
                
                % assigning Electric Conductivity
                obj.edata.ElectricConductivity(obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).ElectricConductivity.value;
                
                % assigning Internal Current Density
                if mzptr.props.isExcited
                    
                    switch mzptr.props.excitation.type
                        
                        case 'currentDensity'
                            
                            obj.edata.InternalCurrentDensity(obj.m.ezi(:, mzptr.zi)) = ...
                                mzptr.props.excitation.value;
                            
                        case 'current'
                            
                            obj.edata.InternalCurrentDensity(obj.m.ezi(:, mzptr.zi)) = ...
                                mzptr.props.excitation.value / mzptr.getArea;
                            
                    end
                    
                end
                
                % assigning Magnetization
                if mzptr.props.isMagnetized
                    
                    M = mzptr.props.magnetization.getM(mzptr.getCenterOfElements);
                    obj.edata.MagnetizationX(obj.m.ezi(:, mzptr.zi)) = M(:, 1)';
                    obj.edata.MagnetizationY(obj.m.ezi(:, mzptr.zi)) = M(:, 2)';
                    
                end
                
            end
            
            % applying current of excitation matrices
            mNames = fieldnames(obj.exmtcs);
            
            for i = 1:numel(mNames)
                
                % matrix pointer
                mptr = obj.exmtcs.(mNames{i});

                for j = 1:mptr.Nmzs
                    
                    mzptr = obj.m.mzs.(mptr.mzsName{j});
                    cptr = obj.coils.(mptr.mzsName{j});
                    obj.edata.InternalCurrentDensity(obj.m.ezi(:, mzptr.zi)) = ...
                        cptr.sign * cptr.turns * mptr.current / mptr.np / mzptr.getArea;
                    
                end
                
            end
            
            disp('Initialization of material and force data compeleted.')
            toc, disp('-------------------------------------------------------');
            
            % change states
            obj.isElementDataAssigned = true;
            
        end
        
        function solve(obj)
            
            % prerequisties
            obj.assignEdata;
            
            % updating boundary conditions
            obj.bcs.updateAll;
            
            % getting mesh zone names
            mzNames = fieldnames(obj.m.mzs);
            
            % Construction of [K] and [F]
            tic, disp('-------------------------------------------------------');
            
            % Assembeling [F]
            fi = repmat(obj.edata.InternalCurrentDensity * obj.units.k_currentDensity, 3, 1);
            Mx = repmat((obj.edata.MagnetizationX .* obj.m.gea) * obj.units.k_magnetisation, 3, 1);
            My = repmat((obj.edata.MagnetizationY .* obj.m.gea) * obj.units.k_magnetisation, 3, 1);
            F = (fi .* repmat((obj.m.gea / 3), 3, 1) + (obj.m.mtcs.gphiy .* Mx - obj.m.mtcs.gphix .* My) / obj.units.k_length);
            
            % applying scales on load vector
            F = F * obj.units.k_length^2 / obj.units.k_magneticVectorPotential;
            F = sparse(obj.m.cl', ones(3 * obj.m.Ne, 1), F);
            
            % Assembeling [K]
            [Iindex, Jindex] = getij(3, 1);
            K = sparse(obj.m.cl(:, Iindex)', obj.m.cl(:, Jindex)', ...
                repmat(obj.edata.MagneticReluctivity, 9, 1) .* ...
                obj.m.mtcs.Ke(getkindex(3), :));
            disp('Construction of [K] and [F] compeleted.');
            toc, disp('-------------------------------------------------------');
            tic, disp('-------------------------------------------------------');
            
            % imposing boundary conditions on [K] and [F]
            % dbcs
            if obj.bcs.Nd
                F(obj.bcs.iD) = obj.bcs.vD;
                K(obj.bcs.iD, :) = sparse(1:obj.bcs.Ndbcs, obj.bcs.iD, ones(1, obj.bcs.Ndbcs), obj.bcs.Ndbcs, obj.m.Nn);
            end
            
            % opbcs
            if obj.bcs.Nop
                F(obj.bcs.mOP) = F(obj.bcs.mOP) - F(obj.bcs.sOP);
                F(obj.bcs.sOP) = 0;
                K(obj.bcs.mOP, :) = K(obj.bcs.mOP, :) - K(obj.bcs.sOP, :);
                K(obj.bcs.sOP, :) = sparse([1:obj.bcs.Nopbcs, 1:obj.bcs.Nopbcs], ...
                    [obj.bcs.mOP; obj.bcs.sOP], ones(1, 2 * obj.bcs.Nopbcs), obj.bcs.Nopbcs, obj.m.Nn);
            end
            
            % epbcs
            if obj.bcs.Nep
                F(obj.bcs.mEP) = F(obj.bcs.mEP) + F(obj.bcs.sEP);
                F(obj.bcs.sEP) = 0;
                K(obj.bcs.mEP, :) = K(obj.bcs.mEP, :) + K(obj.bcs.sEP, :);
                K(obj.bcs.sEP, :) = sparse([1:obj.bcs.Nepbcs, 1:obj.bcs.Nepbcs], ...
                    [obj.bcs.mEP; obj.bcs.sEP], [ones(1, obj.bcs.Nepbcs), -ones(1, obj.bcs.Nepbcs)], obj.bcs.Nepbcs, obj.m.Nn);
            end
            
            disp('All boundary condition imposed.');
            toc, disp('-------------------------------------------------------');
            
            % solving [K][U] = [F]
            tic, disp('-------------------------------------------------------');
            
            % solving equation KU = F
            if ~any(F)
                obj.results.A = full(F);
                return
            end
            
            obj.results.A = full(K \ F);
            obj.evalBe;
            disp('initial geuss evaluated.')
            toc, disp('-------------------------------------------------------');
            
            % loop for nonlinear solver
            tic, disp('-------------------------------------------------------');
            
            % initials values
            RelError = inf;
            Iterations = 0;
            xNgt = obj.m.Ne;
            xNgp = obj.m.Nn;
            
            % preparing error monitoring
            if obj.monitorResiduals

                ERF = gcf; cla;
                set(ERF, 'Name', 'md2d_tl3_ihnl solver', ...
                    'WindowStyle', 'Normal', 'ToolBar', 'none', 'Menu', 'none');
                er = animatedline('color', 'r', 'Linewidth', 1.2, 'Marker', 's', 'MarkerEdgeColor','k');
                title('Relative Error');
                ylabel('log10(|dA|/|A|)');
                xlabel('Iteration Number');
                set(gca, 'box', 'on');

            end
            
            % solver history
            obj.solverHistory.relativeError = [];
            obj.solverHistory.totalEnergy = [];
            obj.solverHistory.totalConergy = [];
            
            alpha = linspace(0.8, 0.95, obj.solverSettings.maxIteration);
            
            % loop for non-linearity
            fprintf('Iter|Error   |Residual|time\n');
            while (RelError > obj.solverSettings.relativeError) && (Iterations < obj.solverSettings.maxIteration)
                
                % starting loop time
                loopTime = tic;
                
                % evaluation of B2 for each elements
                obj.evalBe;
                [obj.solverHistory.totalEnergy(end + 1),obj.solverHistory.totalConergy(end + 1)] = obj.evalTotalEnergyCoenergy;
                Bk = sqrt(obj.results.Bex.^2 + obj.results.Bey.^2);
                
                % updating nu
                for i = 1:obj.m.Nmzs
                    mzptr = obj.m.mzs.(mzNames{i});
                    
                    if ~obj.m.mts.(mzptr.material).MagneticPermeability.isLinear
                        obj.edata.MagneticReluctivity(obj.m.ezi(:, mzptr.zi)) = ...
                            ppval(obj.m.mts.(mzptr.material).vB, ...
                            Bk(obj.m.ezi(:, mzptr.zi)));
                    end
                    
                end
                
                % updating dnudB2
                dnudB2 = zeros(1, xNgt);
                
                for i = 1:obj.m.Nmzs
                    mzptr = obj.m.mzs.(mzNames{i});
                    
                    if ~obj.m.mts.(mzptr.material).MagneticPermeability.isLinear
                        dnudB2(obj.m.ezi(:, mzptr.zi)) = ...
                            ppval(obj.m.mts.(mzptr.material).dvdB, ...
                            Bk(obj.m.ezi(:, mzptr.zi))) ./ ...
                            (2 * Bk(obj.m.ezi(:, mzptr.zi)));
                    end
                    
                end
                
                % construction of stiffness matrix [K]
                K = sparse(obj.m.cl(:, Iindex)', ...
                    obj.m.cl(:, Jindex)', ...
                    repmat(obj.edata.MagneticReluctivity, 9, 1) .* ...
                    obj.m.mtcs.Ke(getkindex(3), :));
                
                % construction of [K] and [F] in NR algorithm
                FF = -K * obj.results.A + F;
                
                % evaluation and adding of jacobian matrix
                K = K + sparse(obj.m.cl(:, Iindex)', obj.m.cl(:, Jindex)', ...
                    IHNLNRMSTL3_evalG(obj.m.cl, obj.m.mtcs.Ke, obj.m.mtcs.gphix, ...
                    obj.m.mtcs.gphiy, obj.results.A, dnudB2) / obj.units.k_length^2);
                
                % imposing boundary conditions on incrimentals
                % dbcs
                if obj.bcs.Nd
                    FF(obj.bcs.iD) = 0;
                    K(obj.bcs.iD, :) = sparse(1:obj.bcs.Ndbcs, obj.bcs.iD, 1, obj.bcs.Ndbcs, xNgp);
                end
                
                % opbcs
                if obj.bcs.Nop
                    FF(obj.bcs.mOP) = FF(obj.bcs.mOP) - FF(obj.bcs.sOP);
                    FF(obj.bcs.sOP) = 0;
                    K(obj.bcs.mOP, :) = K(obj.bcs.mOP, :) - K(obj.bcs.sOP, :);
                    K(obj.bcs.sOP, :) = sparse([1:obj.bcs.Nopbcs, 1:obj.bcs.Nopbcs], ...
                        [obj.bcs.mOP; obj.bcs.sOP], ones(1, 2 * obj.bcs.Nopbcs), obj.bcs.Nopbcs, xNgp);
                end
                
                % epbcs
                if obj.bcs.Nep
                    FF(obj.bcs.mEP) = FF(obj.bcs.mEP) + FF(obj.bcs.sEP);
                    FF(obj.bcs.sEP) = 0;
                    K(obj.bcs.mEP, :) = K(obj.bcs.mEP, :) + K(obj.bcs.sEP, :);
                    K(obj.bcs.sEP, :) = sparse([1:obj.bcs.Nepbcs, 1:obj.bcs.Nepbcs], ...
                        [obj.bcs.mEP; obj.bcs.sEP], [ones(1, obj.bcs.Nepbcs), -ones(1, obj.bcs.Nepbcs)], obj.bcs.Nepbcs, xNgp);
                end
                
                % solving [K][U] = [F]
                dU = K \ FF;
                obj.results.A = full(obj.results.A + alpha(Iterations+1)*dU);
                
                % check for convergency
                Residual = norm(dU, 1);
                RelError = Residual / norm(obj.results.A, 1);
                
                % monitoring of error
                if obj.monitorResiduals
                    addpoints(er, Iterations, log10(RelError));
                    drawnow;
                end
                
                % solver history
                obj.solverHistory.relativeError(end + 1) = RelError;
                
                % printing Residual and RelError
                fprintf('->%2d|%.2e|%.2e|%0.3f\n', Iterations, RelError, Residual, toc(loopTime));
                
                % go to next iteration
                Iterations = Iterations + 1;
            end
            
            obj.solverHistory.iterations = Iterations;
            disp(['Number of total iterations = ', num2str(Iterations - 1)]);
            toc, disp('-------------------------------------------------------');
            
            % change states
            obj.evalBe;
            obj.isBnEvaluated = false;
            
        end
        
        % Solver History
        function plotRelativeError(obj)
            
            figure('Name', 'emdlab -> ms2d_ihnl solver');
            plot(1:obj.solverHistory.iterations, ...
                log10(obj.solverHistory.relativeError), 'color', 'r', 'Linewidth', 1.2, ...
                'marker', 's', 'MarkerEdgeColor', 'b');
            title('Relative Error (|dA|/|A|)')
            ylabel('log10(|dA|/|A|)')
            xlabel('Iteration Number')
            if obj.solverHistory.iterations == 1, return; end
            set(gca, 'xlim', [1, obj.solverHistory.iterations]);
            
        end
        
        function plotEnergy(obj)

            figure('Name', 'emdlab -> ms2d_ihnl solver');
            plot(1:obj.solverHistory.iterations, ...
                obj.solverHistory.totalEnergy, 'color', 'r', 'Linewidth', 1.2, ...
                'marker', 's', 'MarkerEdgeColor', 'b');
            title('Total Energy')
            ylabel('Energy [J]')
            xlabel('Iteration Number')
            if obj.solverHistory.iterations == 1, return; end
            set(gca, 'xlim', [1, obj.solverHistory.iterations]);

        end
        
        function plotEnergyResidual(obj)

            figure('Name', '[EMDLAB] IHNLNRMSTL3 Solver', 'NumberTitle', 'off');
            plot(2:obj.solverHistory.iterations, ...
                diff(obj.solverHistory.totalEnergy), 'color', 'r', 'Linewidth', 1.2, ...
                'marker', 's', 'MarkerEdgeColor', 'b');
            title('Residual of Energy')
            ylabel('Residual [J]')
            xlabel('Iteration Number')
            if obj.solverHistory.iterations == 1, return; end
            set(gca, 'xlim', [1, obj.solverHistory.iterations]);
            
        end
        
    end
    
end
