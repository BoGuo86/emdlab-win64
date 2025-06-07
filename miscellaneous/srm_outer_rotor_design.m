classdef srm_outer_rotor_design < handle
    properties
        %% main parameters
        % rated out put power [watt]
        Pout
        % rated output speed [rpm]
        rpm
        % DC supply voltage [volt]
        Vdc
        % desired efficiency [0-1]
        etad
        % peak of phase current [A]
        Ip
        % air gap mechanical loading [Pa]
        sigmaF
        % ld ratio of motor
        LD
        % saturation flux density in aligned position [Tesla]
        Bsat
        % number of stator tooth
        Ns
        % number of rotor tooth
        Nr
        % number of phases
        m
        % torque per rotor volume [N/m^3]
        TRV
        %% lamination data
        % lamination density
        rhoIron
        % stacking factor
        Kst
        % weight of stator iron [Kg]
        statorWeight
        % weight of rotor iron
        rotorWeight
        %% winding parameters
        % coil resistance [ohm]
        Rcoil
        % phase resistance [ohm]
        Rph
        % mean turn length [mm]
        Lmt
        % number of coil per phase
        Ncph
        % number of coil turns
        Ncoil
        % conductivity of winding wire [S/m]
        sigmaWire
        % winding wire diameter
        Dwire
        % number of parallel path
        np
        % number of parallel path
        a
        % winding fill factor
        Kf
        
        
        
        %% dimensions
        wsyPwsp
        wryPwrp
        alphas
        alphar
        g
        
        Rro
        Rso
        
        %
        D
        
        
        
        
        Ts
        Tr
        hr
        hs
        wry
        wsy
        wsp
        wrp
        betar
        betas
        
        rsh
        Rrm
        Rsm
        Nsplit
        Lst
        As
        Asp
        Arp
        Asy
        Ary
        
        Wcon
        
        
        
        
        
    end
    properties (Hidden = true)
        gammar;
        gammas;
        mesh_thetar1 = 8;
        mesh_thetar2 = 4;
        mesh_thetag = 2;
        mesh_thetas1 = 3;
        mesh_thetas2 = 4;
        geom = TMDBC;
    end
    properties (Access = private)
        isSetGeom = false;
    end
    properties (Dependent = true)
        %% mechanical parameters
        % mechanical angular velocity [rad/sec]
        wm
        % ratef output torque [N.m]
        Tout
        % rated rotor speed [rps]
        rps
        %% winding parameters
        % number of series coils in a phase
        nser
        % stroke angle
        e
    end
    methods
        function obj = srm_outer_rotor_design()
        end
        %% Setters
        function set.Pout(obj,value)
            if value<0
                warning('Output power must be a real positive number');
            else
                obj.Pout = value;
            end
        end
        function set.rpm(obj,value)
            if value<0
            else
                obj.rpm = value;
            end
        end
        function set.Rro(obj,value)
            if value<0
            else
                obj.Rro = value;
            end
        end
        function set.Vdc(obj,value)
            if value<0
            else
                obj.Vdc = value;
            end
        end
        function set.Ns(obj,value)
            if value<0 || rem(value,2)
                warning('number of stator tooth must be a positive even integer');
            else
                obj.Ns = value;
            end
        end
        function set.Nr(obj,value)
            if value<0 || rem(value,2)
                warning('number of rotor tooth must be a positive even integer');
            else
                obj.Nr = value;
            end
        end
        function set.Ncoil(obj,value)
            if value<0 || rem(value,1)
                warning('number of coil turns must be a positive integer');
            else
                obj.Ncoil = value;
            end
        end
        function set.m(obj,value)
            if value<0 || rem(value,1)
                warning('number of phases must be a positive integer');
            else
                obj.m = value;
            end
        end
        function set.a(obj,value)
            if value<0 || rem(value,1)
                warning('number of parallel path must be a positive integer');
            else
                obj.a = value;
            end
        end
        function set.D(obj,value)
            if value<0.45 || value>0.85
                warning('Out of range');
            else
                obj.D = value;
            end
        end
        function set.sigmaF(obj,value)
            if value<20e3 || value>40e3
                warning('Out of range');
            else
                obj.sigmaF = value;
            end
        end
        function set.etad(obj,value)
            if value<0.4 || value>1
                warning('Out of range');
            else
                obj.etad = value;
            end
        end
        function set.g(obj,value)
            if value<0.25 || value>0.8
                warning('Out of range');
            else
                obj.g = value;
            end
        end
        function set.Bsat(obj,value)
            if value<1.8 || value>2.4
                warning('Out of range');
            else
                obj.Bsat = value;
            end
        end
        function set.LD(obj,value)
            if value<0.5 || value>3
                warning('Out of range');
            else
                obj.LD = value;
            end
        end
        function set.alphas(obj,value)
            if value<0.2 || value>0.8
                warning('Out of range');
            else
                obj.alphas = value;
            end
        end
        function set.alphar(obj,value)
            if value<0.2 || value>0.8
                warning('Out of range');
            else
                obj.alphar = value;
            end
        end
        function set.wryPwrp(obj,value)
            if value<0.5 || value>1
                warning('Out of range');
            else
                obj.wryPwrp = value;
            end
        end
        function set.wsyPwsp(obj,value)
            if value<0.5 || value>1
                warning('Out of range');
            else
                obj.wsyPwsp = value;
            end
        end
        function set.Kst(obj,value)
            if value<0.9 || value>1
                warning('Out of range');
            else
                obj.Kst = value;
            end
        end
        function set.Kf(obj,value)
            if value<0.2 || value>0.8
                warning('Out of range');
            else
                obj.Kf = value;
            end
        end
        %% Getters: dependent variables
        function y = get.Ts(obj)
            y = 2*pi/obj.Ns;
        end
        function y = get.Tr(obj)
            y = 2*pi/obj.Nr;
        end
        function y = get.Ip(obj)
            y = obj.Pout/obj.etad/obj.Vdc;
        end
        function y = get.Nsplit(obj)
            y = gcd(obj.Ns,obj.Nr);
        end
        function y = get.Lmt(obj)
            y = obj.Lst+2*pi*obj.wsp;
        end
        function y = get.Rcoil(obj)
            y = obj.Ncoil^2*obj.Lmt/obj.sigmaWire/obj.As/2/obj.Kf/1e-3;
        end
        function y = get.Ncph(obj)
            y = obj.Ns/obj.m;
        end
        function y = get.nser(obj)
            y = obj.Ncph/obj.np;
        end
        function y = get.e(obj)
            y = 2*pi/obj.m/obj.Nr;
        end
        function y = get.rps(obj)
            y = obj.rpm/60;
        end
        function y = get.wm(obj)
            y = 2*pi*obj.rps;
        end
        function y = get.Tout(obj)
            y = obj.Pout/obj.wm;
        end
        function y = get.Rph(obj)
            y = obj.Rcoil*obj.nser/obj.np;
        end
        function writeParFile(obj)
            writeParFile('Rso',obj.Rso,'Rro',obj.Rro,'hs',obj.hs,...
                'wsy',obj.wsy,'g',obj.g,'Ts',obj.Ts,'Tr',obj.Tr,...
                'wry',obj.wry,'hr',obj.hr,'gammas',obj.gammas,...
                'gammar',obj.gammar,'mesh_thetag',obj.mesh_thetag,...
                'mesh_thetar1',obj.mesh_thetar1,...
                'mesh_thetar2',obj.mesh_thetar2,...
                'mesh_thetas1',obj.mesh_thetas1,...
                'mesh_thetas2',obj.mesh_thetas2,...
                'betas',obj.betas,'betar',obj.betar,...
                'rsh',obj.rsh);
        end
        %% evaluations
        % Lst, Rso, LD and sigmaF
        function eval_Lst_Rso_For_sigmaF_LD(obj)
            tmp = obj.Tout*2/(pi*obj.sigmaF*obj.LD);
            obj.Rso = nthroot(tmp,3)*1e3 - obj.g;
            obj.Rso = round(obj.Rso,1);
            obj.Lst = 2*(obj.Rso+obj.g)*obj.LD;
            obj.Lst = round(obj.Lst,1);
        end
        %% updates
        % Lst, Rso and LD
        function update_Lst_For_LD_Rso(obj)
            obj.Lst = obj.LD*2*(obj.Rso+obj.g);
        end
        function update_LD_For_Lst_Rso(obj)
            obj.LD = obj.Lst/(2*(obj.Rso+obj.g));
        end
        function update_Rso_For_LD_Lst(obj)
            obj.Rso = obj.Lst/obj.LD/2 - obj.g;
        end
        % D, Rro and Rso
        function update_D_For_Rso_Rro(obj)
            obj.D = obj.Rso/obj.Rro;
        end
        function update_Rro_For_D_Rso(obj)
            obj.Rro = obj.Rso/obj.D;
        end
        function update_Rso_For_D_Rro(obj)
            obj.Rso = obj.Rro*obj.D;
        end
        % betas and alphas
        function update_betas_For_alphas(obj)
            obj.betas = obj.alphas*obj.Ts;
        end
        function update_alphas_For_betas(obj)
            obj.alphas = obj.betas/obj.Ts;
        end
        % betar and alphar
        function update_betar_For_alphar(obj)
            obj.betar = obj.alphar*obj.Tr;
        end
        function update_alphar_For_betar(obj)
            obj.alphar = obj.betar/obj.Tr;
        end
        % wsp and betas
        function update_wsp_For_betas(obj)
            obj.wsp = 2*(obj.Rso)*sin(obj.betas);
        end
        function update_betas_For_wsp(obj)
            obj.betas = asin(obj.wsp/2/obj.Rso);
        end
        % wrp and betar
        function update_wrp_For_betar(obj)
            obj.wrp = 2*(obj.Rso+obj.g)*sin(obj.betar);
        end
        function update_betar_For_wrp(obj)
            obj.betar = asin(obj.wrp/2/(obj.Rso+obj.g));
        end
        % wsy and wsyPwsp
        function update_wsy_For_wsyPwsp(obj)
            obj.wsy = obj.wsyPwsp*obj.wsp;
        end
        function update_wsyPwsp_For_wsy(obj)
            obj.wsyPwsp = obj.wsy/obj.wsp;
        end
        % wry and wryPwrp
        function update_wry_For_wryPwrp(obj)
            obj.wry = obj.wryPwrp*obj.wrp;
        end
        function update_wryPwrp_For_wry(obj)
            obj.wryPwrp = obj.wry/obj.wrp;
        end
        %% function in flowchart
        function evalMainDimensions(obj)
            obj.eval_Lst_Rso_For_sigmaF_LD;
            obj.update_Rro_For_D_Rso;
            obj.update_betar_For_alphar;
            obj.update_betas_For_alphas;
            obj.update_wrp_For_betar;
            obj.update_wsp_For_betas;
            obj.update_wsy_For_wsyPwsp;
            obj.update_wry_For_wryPwrp;
        end
        %% check constraints
        function y = checkPoleAngles(obj)
            if min(obj.betar,obj.betas) <= obj.e
                warning('Stroke angle must be grater than betar and betas.');
                y = false;
            end
        end
    end
end