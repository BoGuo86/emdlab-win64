classdef ID_SRM0<handle
    properties (SetAccess = private)
        Pout
        rpm
        Vdc 
        Rso 
        D 
        g 
        wsyPwsp
        wryPwrp 
        hrPg 
        alphas
        alphar 
        ns
        nr    
        Bsat 
        win = struct();
        LD 
        np 
        Kf 
        sigmaWin = 58e6; 
        rhoIron = 7850; 
        Kst 
        mesh_thetar1 = 5;
        mesh_thetar2 = 4;
        mesh_thetag = 2;
        mesh_thetas1 = 3;
        mesh_thetas2 = 4;
    end
    properties (Dependent=true)
        Rro
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
        gammar
        gammas
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
        nser
        Wiron
        Rcoil
        Rph
        Lmt
        Ncph
        rps
        wm
        T
        e
        Wcon
        Np
    end
    methods
        function obj = ID_SRM0()
            obj.set_Pout(500);
            obj.set_rpm(400);
            obj.set_ns(8);
            obj.set_nr(6);
            obj.set_Vdc(310);
            obj.set_Rso(50);
            obj.set_D(0.6);
            obj.set_alphas(0.5);
            obj.set_alphar(0.5);
            obj.set_g(0.35);
            obj.set_LD(0.8);
            obj.set_hrPg(25);
            obj.set_wsyPwsp(0.6);
            obj.set_wryPwrp(0.6);
            obj.set_Kst(0.97);
            obj.set_Kf(0.45);
            obj.set_np(1);
            obj.set_Bsat(2.1);
        end
        function set_Pout(obj,value)
            if value<0
            else
                obj.Pout = value;
            end          
        end
        function set_rpm(obj,value)
            if value<0
            else
                obj.rpm = value;
            end          
        end
        function set_Rso(obj,value)
            if value<0
            else
                obj.Rso = value;
            end          
        end
        function set_Vdc(obj,value)
            if value<0
            else
                obj.Vdc = value;
            end          
        end
        function set_ns(obj,value)
            if value<0 || rem(value,2)
            else
                obj.ns = value;
            end          
        end
        function set_nr(obj,value)
            if value<0 || rem(value,2)
            else
                obj.nr = value;
            end          
        end
        function set_np(obj,value)
            if value<0 || rem(value,1)
            else
                obj.np = value;
            end          
        end
        function set_D(obj,value)
            if value<0.45 || value>0.85
            else
                obj.D = value;
            end          
        end
        function set_g(obj,value)
            if value<0.25 || value>0.8
            else
                obj.g = value;
            end          
        end
        function set_Bsat(obj,value)
            if value<1.8 || value>2.4
            else
                obj.Bsat = value;
            end          
        end
        function set_LD(obj,value)
            if value<0.5 || value>3
            else
                obj.LD = value;
            end          
        end
        function set_alphas(obj,value)
            if value<0.2 || value>0.8
            else
                obj.alphas = value;
            end          
        end
        function set_alphar(obj,value)
            if value<0.2 || value>0.8
            else
                obj.alphar = value;
            end          
        end
        function set_wryPwrp(obj,value)
            if value<0.5 || value>1
            else
                obj.wryPwrp = value;
            end          
        end
        function set_wsyPwsp(obj,value)
            if value<0.5 || value>1
            else
                obj.wsyPwsp = value;
            end          
        end
        function set_Kst(obj,value)
            if value<0.9 || value>1
            else
                obj.Kst = value;
            end          
        end
        function set_Kf(obj,value)
            if value<0.2 || value>0.8
            else
                obj.Kf = value;
            end          
        end
        function set_hrPg(obj,value)
            if value<20 || value>30
            else
                obj.hrPg = value;
            end          
        end
        function y = get.Rro(obj)
            y = obj.Rso*obj.D;
        end
        function y = get.Ts(obj)
            y = 2*pi/obj.ns;
        end
        function y = get.Tr(obj)
            y = 2*pi/obj.nr;
        end
        function y = get.hr(obj)
            y = obj.hrPg*obj.g;
        end
        function y = get.betas(obj)
            y = obj.alphas*obj.Ts;
        end
        function y = get.betar(obj)
            y = obj.alphar*obj.Tr;
        end
        function y = get.wsp(obj)
            y = 2*(obj.Rro+obj.g)*sin(obj.betas);
        end
        function y = get.wrp(obj)
            y = 2*obj.Rro*sin(obj.betar);
        end
        function y = get.wsy(obj)
            y = obj.wsyPwsp*obj.wsp;
        end
        function y = get.wry(obj)
            y = obj.wryPwrp*obj.wrp;
        end
        function y = get.rsh(obj)
            y = obj.Rro-obj.hr-obj.wry;
        end
        function y = get.gammas(obj)
            y = asin((obj.Rro+obj.g)*sin(obj.betas/2)/(obj.Rro+obj.g+obj.hs));
        end
        function y = get.gammar(obj)
            y = asin(obj.Rro*sin(obj.betar/2)/(obj.Rro-obj.hr));
        end
        function y = get.hs(obj)
            y = obj.Rso-obj.Rro-obj.g-obj.wsy;
        end
        function y = get.Rrm(obj)
            y = obj.Rro-obj.hr;
        end
        function y = get.Rsm(obj)
            y = obj.Rso-obj.wsy;
        end
        function y = get.Nsplit(obj)
            y = gcd(obj.ns,obj.nr);
        end
        function obj = readwin(obj,FileDir)
            % open and read file
            f = fopen(FileDir,'r'); 
            if ~strcmpi(rmspaces(fgetl(f)),'srm0')
                error('Motor type is wrong ...');
            end
            Data = fscanf(f,'%d');
            fclose(f);
            obj.win.m = Data(2);
            Npp = obj.ns/obj.win.m;
            obj.np = Data(3);
            Data = Data(4:end);
            Data = reshape(Data,2,[])';

            tmp = zeros(Npp,2);
            
            for i = 1:obj.win.m
                for j = 1:Npp
                    tmp(j,1) = Data((i-1)*Npp+j,1);
                    tmp(j,2) = Data((i-1)*Npp+j,2);
                end
                obj.win.(['phase',num2str(i)]) = tmp; 
            end   
        end
        function y = get.Lst(obj)
            y = 2*obj.Rro*obj.LD;
        end
        function y = get.Asp(obj)
            y = obj.wsp*obj.hs;
        end
        function y = get.Arp(obj)
            y = obj.wrp*obj.hr;
        end
        function y = get.As(obj)
            y = pi*(obj.Rsm^2-(obj.Rro+obj.g)^2);
            y = y-obj.ns*obj.Asp;
        end
        function y = get.Asy(obj)
            y = pi*(obj.Rso^2-obj.Rsm^2);
        end
        function y = get.Ary(obj)
            y = pi*(obj.Rrm^2-obj.rsh^2);
        end
        function y = get.Wiron(obj)
            y = obj.Lst*(obj.Ary+obj.Arp+obj.Asp+obj.Asy)*1e-9*obj.rhoIron;
        end
        function y = get.Lmt(obj)
            y = obj.Lst+2*pi*obj.wsp;
        end
        function y = get.Rcoil(obj)
            y = obj.Np^2*obj.Lmt/obj.sigmaWin/obj.As/2/obj.Kf/1e-3; 
        end
        function y = get.Ncph(obj)
            y = obj.ns/obj.win.m;
        end
        function y = get.nser(obj)
            y = obj.Ncph/obj.np;
        end
        function y = get.e(obj)
            y = 2*pi/obj.win.m/obj.nr;
        end
        function y = get.rps(obj)
            y = obj.rpm/60;
        end
        function y = get.wm(obj)
            y = 2*pi*obj.rps;
        end
        function y = get.T(obj)
            y = obj.Pout/obj.wm;
        end
        function y = get.Wcon(obj)
            y = obj.Pout*2*pi/obj.win.m/obj.nr/obj.wm;
        end
        function y = get.Np(obj)
            y = obj.Vdc*obj.betas/(obj.wm*obj.Bsat*obj.Kst*obj.wsp*...
                obj.Lst*1e-6*obj.nser);
            y = ceil(y)*4;
        end
        function y = get.Rph(obj)
            y = obj.Rcoil*obj.nser/obj.np;
        end 
    end    
end 