classdef ID_PMSM0
    properties
        % rated output power
        Pout = 500;
        % rated speed
        rpm = 1500;
        % rated voltage of motor
        Vrms = 220/sqrt(3);
        % mean magnetic flux density in air gap
        Bg = 0.87;
        % stator tooth magnetic flux
        Bt = 1.4;
        % stator yoke magnetic flux
        Bsy = 1.4;
        % rotor yoke magnetic flux
        Bry = 1.4;
        % number of stator slots
        Ns = 27;
        % number of magnets
        Nm = 12;
        % magnet height
        hm = 4;
        % air gap length
        g = 1;
        % stack length
        Lst = 100;
        % outer stator radius
        Rso = 50;
        % ratio of Rro to Rso
        D = 0.6;
        % stacking factor
        Kst = 0.97;
        % fill factor
        Kf = 0.6;
        % slot dimensions
        alphas = 0.15;
        hss1 = 1;
        hss2 = 0.5;
        % number of parallel path
        a = 1;
        % number of phases
        m = 3;
        % ratio of magnet arc to pole pitch arc
        betam = 0.9;
        % power factor
        pf = 0.95;
        % desired efficiency
        eff = 0.85;
        % electrical conductivity
        sigma = 57e6;
    end
    properties(Dependent = true)
        % angular mechanical speed
        wm
        % electrical mechanical speed
        we
        % rated output torque
        T
        % outer rotor radius
        Rro
        % shaft radius
        rsh
        % width of stator yoke
        wsy
        % width of rotor yoke
        wry
        % width of stator tooth base
        wtb
        % stator slot pitch
        Ts
        % rotor pole pitch
        Tr
        % Nsplit
        Nsplit
        % slot dimensions
        Tso1
        Tso2
        hss3
        % number of coil turns
        Ncoil
        % number of phase turns
        Nph
        % rated frequency of motor
        f
        % rated phase current
        Irms
        % total air gap flux
        phi_tot
        % pole air gap flux
        phi_p
        % slot area
        As
        % stator yoke area
        Asy
        % rotor yoke area
        Ary
        % teeth area
        At
        % volume of bared wire
        Vwb
    end  
    methods
        function obj = ID_PMSM0()
        end
        function y = get.Nsplit(obj)
            y = gcd(obj.Ns/obj.m,obj.Nm);
        end
        function y = get.wm(obj)
            y = obj.rpm*pi/30;
        end   
        function y = get.wtb(obj)
            y = 2*pi*obj.Rro*obj.Bg/(obj.Ns*obj.Kst*obj.Bt);
            y = round(y,1);
        end
        function y = get.wsy(obj)
            y = pi*obj.Rro*obj.Bg/(obj.Nm*obj.Kst*obj.Bsy);
            y = round(y,1);
        end
        function y = get.wry(obj)
            y = pi*obj.Rro*obj.Bg/(obj.Nm*obj.Kst*obj.Bry);
            y = round(y,1);
        end
        function y = get.we(obj)
            y = obj.wm*obj.Nm/2;
        end   
        function y = get.T(obj)
            y = obj.Pout/obj.wm;
        end   
        function y = get.rsh(obj)
            y = obj.Rro - obj.hm - obj.wry;
        end
        function y = get.Ts(obj)
            y = 2*pi/obj.Ns;
        end
        function y = get.Tr(obj)
            y = 2*pi/obj.Nm;
        end
        function y = get.Tso1(obj)
            y = obj.alphas*obj.Ts;
        end
        function y = get.Tso2(obj)
            y = obj.Ts - 2*asin((obj.wtb/2)/(obj.Rso-obj.wsy));
        end
        function y = get.hss3(obj)
            y = obj.Rso - obj.wsy - obj.Rro - obj.g - obj.hss1 - ...
                obj.hss1*cos((obj.Ts-obj.Tso1)/2);
        end    
        function y = get.f(obj)
            y = obj.rpm*obj.Nm/120;
        end      
        function y = getRph(obj,As,Span)
            Aw = As/2;
            Rs = obj.Ncoil^2*obj.Lst/Aw/obj.Kf/obj.sigma/1e-3;
            Rave = 0.5*(obj.Rso+obj.Rso);
            L2 = Rave*Span*obj.Ts;
            Lew = L2+4;
            Rew = obj.Ncoil^2*Lew/Aw/obj.Kf/obj.sigma/1e-3;
            Rcoil = 2*(Rs+Rew);
            y = obj.Ns*Rcoil/obj.a^2/obj.m;
        end 
        function y = get.Rro(obj)
            y = obj.D*obj.Rso;
        end
        function y = get.phi_tot(obj)
            y = 2*pi*obj.Rro*obj.Bg*obj.Lst*1e-6;
        end
        function y = get.phi_p(obj)
            y = obj.phi_tot/obj.Nm;
        end
        function y = get.Ncoil(obj)
            y = obj.Vrms/(4.44*obj.f*obj.phi_p*0.955);
            y = obj.a*obj.m*y/obj.Ns;
            y = ceil(y);
        end
        function y = get.Nph(obj)
            y = obj.Ns*obj.Ncoil/(obj.a*obj.m);
        end
        function y = get.Irms(obj)
            Pin = obj.Pout/obj.eff;
            y = Pin/obj.m/obj.Vrms/obj.pf;
        end
        function y = get.Asy(obj)
            y = pi*(obj.Rso^2-(obj.Rso-obj.wsy)^2);
        end
        function y = get.At(obj)
            y = obj.Ns*obj.wtb*(obj.hss3);
        end
        function y = get.As(obj)
            y = pi*((obj.Rso-obj.wsy)^2-(obj.Rso-obj.wsy-obj.hss3)^2);
            y = (y-obj.At)/obj.Ns;
        end
        function y = get.Vwb(obj)
            y = obj.As*obj.Kf*obj.Lst/obj.Ncoil;
        end
        function y = getD(obj)    
            xD = 0.4:0.001:0.7;
            Km = zeros(1,length(xD));
            for i = 1:length(xD)
                obj.D = xD(i);
                Km(i) = obj.Bg*obj.Rro*sqrt(obj.Vwb)/sqrt(obj.sigma);
            end
            plot(xD,Km)
            [~,y] = max(Km);
            y = xD(y);            
        end
    end
end