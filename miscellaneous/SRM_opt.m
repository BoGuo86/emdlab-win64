function SRM_opt()

% creation of result file
f = fopen('Results.txt','w');
fprintf(f, 'Pole shape optimization\n');
fprintf(f, '============================================================\n');
Iterations = 0;
fmincon(@fem_routine,[40,22.5,20],[0,-1,1],0,[],[],[40,15,15],[55,30,30],[],optimset('MaxIter',200));
fprintf(f, '============================================================\n');
fclose(f);

    function Objective = fem_routine(x)
        %% initialization
        materialdir = [cd,'\MaterialsData'];
        %% Mesh
        m = TMDBC;
        m.addMaterial(materialdir,'air');
        m.addMaterial(materialdir,'m19_24ga');
        m.addMaterial(materialdir,'copper');
        % Rotor
        Nt = 6;
        glib_srm_rt1(15,30.3,x(1)-0.36/2,Nt,x(2)*pi/180,15,6,0.45);
        m.read_g2d_bin('geom.g2d','MG1',10);
        m.setMaterial('r1','m19_24ga');
        m.setmzc('r1',[0.7883,0.8529,0.4856]);
        m.setmzc('rap11',[0.9700,0.6951,0.2733]);
        m.cmmz('r2','r1',[1,0]);
        m.cmmz('rap21','rap11',[1,0]);
        for i = 1:2:2*(Nt-1)
            m.crmz(['r',num2str(i+2)],['r',num2str(i)],2*pi/Nt);
            m.crmz(['r',num2str(i+3)],['r',num2str(i+1)],2*pi/Nt);
        end
        for i = 1:Nt-1
            m.crmz(['rap1',num2str(i+1)],['rap1',num2str(i)],2*pi/Nt);
            m.crmz(['rap2',num2str(i+1)],['rap2',num2str(i)],2*pi/Nt);
        end
        temp = getlist('r',1:2*Nt);
        m.jmzs('rotor',temp{:});
        % Stator and Coil
        Nt = 8;
        glib_srm_st1(x(1)+0.36/2,78.4,89.8,Nt,x(3)*pi/180,0.45,2,3);
        m.read_g2d_bin('geom.g2d','MG1',10);
        m.setMaterial('s1','m19_24ga');
        m.setMaterial('c11','copper');
        m.setmzc('s1',[0.8758,0.8650,0.3552]);
        m.cmmz('s2','s1',[1,0]);
        m.cmmz('c21','c11',[1,0]);
        m.setmzc('c11',[0.8921,0.4059,0.6044]);
        m.setmzc('c21','c');
        for i = 1:2:2*(Nt-1)
            m.crmz(['s',num2str(i+2)],['s',num2str(i)],2*pi/Nt);
            m.crmz(['s',num2str(i+3)],['s',num2str(i+1)],2*pi/Nt);
        end
        for i = 1:Nt-1
            m.crmz(['c1',num2str(i+1)],['c1',num2str(i)],2*pi/Nt);
            m.crmz(['c2',num2str(i+1)],['c2',num2str(i)],2*pi/Nt);
        end
        temp = getlist('s',1:2*Nt);
        m.jmzs('stator',temp{:});
        m.ggmesh;
        kr = m.getnIndexOnCircle([0,0],x(1)-0.36/2);
        ks = m.getnIndexOnCircle([0,0],x(1)+0.36/2);
        rps = m.nodes(kr,:);
        sps = m.nodes(ks,:);
        %% Solver
        s = IHNLNRMSTL3(m);clear m
        s.setUnit('length', 'mm');
        s.setUnit('currentDensity', 'A/mm^2');
        s.setUnit('magneticVectorPotential', 'A/m');
        % coils of phase a
        Ncoil = 56;
        s.setDepth(151);
        pos_a = {'c11','c25'};
        neg_a = {'c21','c15'};
        %% proccess
        phaseCurrent = linspace(0,20,2);
        rotorAngle = linspace(0,pi/6,20);
        Ni = length(phaseCurrent);
        Ntheta = length(rotorAngle);
        linkageFlux = zeros(Ntheta,Ni);
        theta_old = 0;
        Arriko = zeros(Ntheta,Ni);
        % loop for calculation of linkage flux
        for i = 1:Ntheta
            % rotation of regions
            s.m.rmz('rotor',rotorAngle(i)-theta_old);
            for j = 1:6
                s.m.rmz(['rap1',num2str(j)],rotorAngle(i)-theta_old);
                s.m.rmz(['rap2',num2str(j)],rotorAngle(i)-theta_old);
            end
            s.addmz('AG',getmz(rotaterps(MC_AG(sps,rps),rotorAngle(i)),4));
            theta_old = rotorAngle(i);
            s.m.ggmesh;
            % getting index fo boundary conditions
            s.bcs.clearAllBCs;
            s.bcs.setDirichlet(s.m.getfbn,0);
            for j = 2:Ni
                % set excitation
                for k = 1:2
                    s.setExcitation(pos_a{k},2.5, 'cd');
                    s.setExcitation(neg_a{k},-2.5, 'cd');
                end
                s.solve;
                % eval linkage flux
                for k = 1:2
                    linkageFlux(i,j) = linkageFlux(i,j) + s.evallf(pos_a{k}) - s.evallf(neg_a{k});
                end
                Arriko(i,j) = s.evalTorqueByArriko(0.36);
            end
            s.removemz('AG');
        end
        delete(s);
        
        Objective = mean(Arriko(:,2));
        plot(Arriko(:,2)); 
        Iterations = Iterations + 1;
        fprintf(f, '%2d:\tRAIRGAP=%f\tbetar=%f\tbetas=%f\tTorque=%f\n', Iterations, x(1), x(2), x(3), Objective);
        Objective = -Objective;
        
    end

end
