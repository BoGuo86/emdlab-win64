classdef SD_PMSM0
    properties
        % rated output power
        Pout = 500;
        % rated speed
        rpm = 1500;
        % rated voltage of motor
        Vrms = 220;
        % number of magnets
        Nm = 15;
        % number of stator slots
        Ns = 4;
        % magnet height
        hm = 4;
        % air gap length
        g = 1;
        % stack length
        Lst = 100;
        % outer stator radius
        Rso = 50;
        % outer rotor radius
        Rro = 30;
        % stacking factor
        Kst = 0.97;
        % fill factor
        Kf = 0.45;
        % slot dimensions
        bss1 = 0.3;
        hss1 = 1.5;
        hss2 = 0.5;
        % number of parallel path
        a = 1;
        % number of phases
        m = 3;
        % ratio of magnet arc to pole pitch arc
        betam = 0.8;
        % power factor
        pf = 0.95;
        % desired efficiency
        eff = 0.85;
        % width of stator yoke
        wsy
        % width of rotor yoke
        wry
        % width of stator tooth base
        wtb
        % number of coil turns
        Ncoil
        % electrical conductivity
        sigma = 57e6;
        % average winding temperature
        winT = 60;
        % alphaT
        alphaT = 0.1;
        dwb
        nd
    end
    properties(Dependent = true)
        % angular mechanical speed
        wm
        % electrical mechanical speed
        we
        % rated output torque
        T
        % shaft radius
        rsh
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
        % number of phase turns
        Nph
        % rated frequency of motor
        f
        % rated phase current
        Irms
        % volume of stator yoke
        Vsy
        % volume of stator tooth
        Vst
        % area of stator slot
        As
    end  
    methods
        function obj = SD_PMSM0()
        end
        function y = get.Nsplit(obj)
            y = gcd(obj.Ns/obj.m,obj.Nm);
        end
        function y = get.wm(obj)
            y = obj.rpm*pi/30;
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
            y = asin(obj.bss1/2/(obj.Rro+obj.g));
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
        function y = get.Nph(obj)
            y = obj.Ns*obj.Ncoil/(obj.a*obj.m);
        end
        function y = get.Vsy(obj)
            r1 = obj.Rso - obj.wsy;
            r2 = obj.Rso;
            y = pi*(r2^2-r1^2)*obj.Lst;
        end
        function y = get.Vst(obj)
            y = obj.Ns*obj.wtb*obj.hss3*obj.Lst;
        end
        
        function y = getKR(obj)
            
            delta = 1/sqrt(pi*obj.f*obj.sigma*4*pi*1e-7);
            
            Delta = obj.dwb*1e-3/delta;
            
            F = (sinh(2*Delta)+sin(2*Delta))/(cosh(2*Delta)-cos(2*Delta));
            G = (sinh(Delta)-sin(Delta))/(cosh(Delta)+cos(Delta));
            
            y = Delta*F+(2/3)*Delta*G*(obj.nd^2-1);
            
        end
        function readFromFile(obj,Dir)
            File = fopen(Dir);
            if ~File
                error('Wrong file directory');
            end
 
            while true
                str = fgets(File);
                if str==-1
                    break
                end
                str = strsplit(str);
                if isfield(str{1},obj)
                    obj.str{1} = str2double(str{2});
                end
            end
        end
    end
end