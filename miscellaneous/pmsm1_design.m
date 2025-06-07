classdef pmsm1_design < handle
    properties (SetAccess = private)
        % rated output power
        Pout
        % rated speed
        rpm
        % rated voltage of motor
        Vrms
        % mean magnetic flux density in air gap
        Bg
        % stator tooth magnetic flux
        Bt
        % stator yoke magnetic flux
        Bsy
        % rotor yoke magnetic flux
        Bry
        % number of stator slots
        Ns
        % number of magnets
        Nm
        % magnet height
        hm
        % air gap length
        g
        % stack length
        Lst
        % outer stator radius
        Rso
        % ratio of Rso to Rro 
        D
        % stacking factor
        Kst
        % fill factor
        Kf
        % slot dimensions
        bss1
        hss1
        % number of parallel path
        a
        % number of phases
        m
        % power factor
        pf
        % desired efficiency
        eff
        % electrical conductivity
        sigma = 58e6;
        % windinga data
        win = struct();
        % data for mesh
        thetar
        thetag
        thetas    
        % width of stator yoke
        wsy
        % width of rotor yoke
        wry
        % width of stator tooth base
        wtb
        % number of coil turns
        Ncoil
        % Coil Span
        Span
        % end winding length: from stator to bending
        Lend
        winT
        alphaT
        % motor geometry
        geom = TMDBC;
        % Slot area
        Aslot
        % shaft radius
        rsh
        % angular mechanical speed
        wm
        % electrical mechanical speed
        we
        % rated output torque
        T
        % outer rotor radius
        Rro
        % stator slot pitch
        Ts
        % rotor pole pitch
        Tr
        % slot dimensions
        Tso1
        Tso2
        Tso3
        rss1
        rss2
        % number of phase turns
        Nph
        % rated frequency of motor
        f
        % rated phase current
        Irms
        % coil resistance
        Rdc_coil
        Rdc_coil_T
        Rac_coil_T
        kR
        % phase resistance
        Rph
        % CopperLoss
        CopperLoss
        % total air gap flux
        phi_tot
        % pole air gap flux
        phi_p  
    end
    properties (Access = private)
        isSetWin = false;
        isSetGeom = false;
    end
    methods
        function obj = pmsm1_design()
            obj.set_Pout(0.5);
            obj.set_rpm(1500);
            obj.set_Vrms(220/sqrt(3));
            obj.set_Rso(50);
            obj.set_D(0.62);
            obj.set_m(3);
            obj.set_Ns(18);
            obj.set_Nm(10);
            obj.set_Lst(100);
            obj.set_g(1);
            obj.set_hm(4);
            obj.set_betam(0.9);
            obj.set_Kst(0.97);
            obj.set_wry(6);
            obj.set_wsy(6);
            obj.set_wtb(4);
            obj.set_bss1(5);
            obj.set_hss1(2);
            obj.set_Ncoil(17);
            obj.set_Kf(0.4);
            obj.set_eff(0.85);
            obj.set_pf(0.95);
            obj.set_a(1);
            obj.set_thetar(5);
            obj.set_thetag(2);
            obj.set_thetas(4);
            obj.set_Span(3);
            obj.set_Lend(10);
            obj.set_Bg(0.87);
            obj.set_Bt(1.7);
            obj.set_Bry(1.4);
            obj.set_Bsy(1.4);
            obj.set_winT(80);
            obj.set_alphaT(3.8e-3);
            obj.set_kR(1);
        end
        %% Setters
        function BOOL = setVariable(obj,varname,value,varargin)
            switch strrep(varname,' ' , '')
                case 'Pout'
                    BOOL = set_Pout(obj,value,varargin{:});
                case 'rpm'
                    BOOL = set_rpm(obj,value,varargin{:});
                case 'Vrms'
                    BOOL = set_Vrms(obj,value,varargin{:});
                case 'Ns'
                    BOOL = set_Ns(obj,value,varargin{:});
                case 'Nm'
                    BOOL = set_Nm(obj,value,varargin{:});
                case 'g'
                    BOOL = set_g(obj,value,varargin{:});
                case 'Rso'
                    BOOL = set_Rso(obj,value,varargin{:});
                case 'hm'
                    BOOL = set_hm(obj,value,varargin{:});
                case 'Ncoil'
                    BOOL = set_Ncoil(obj,value,varargin{:});
                case 'betam'
                    BOOL = set_betam(obj,value,varargin{:});
                case 'Lst'
                    BOOL = set_Lst(obj,value,varargin{:});
                case 'Kst'
                    BOOL = set_Kst(obj,value,varargin{:});
                case 'm'
                    BOOL = set_m(obj,value,varargin{:});
                case 'D'
                    BOOL = set_D(obj,value,varargin{:});
                case 'wry'
                    BOOL = set_wry(obj,value,varargin{:});
                case 'wsy'
                    BOOL = set_wsy(obj,value,varargin{:});
                case 'wtb'
                    BOOL = set_wtb(obj,value,varargin{:});
                case 'bss1'
                    BOOL = set_bss1(obj,value,varargin{:});
                case 'hss1'
                    BOOL = set_hss1(obj,value,varargin{:});
                case 'Kf'
                    BOOL = set_Kf(obj,value,varargin{:});
                case 'a'
                    BOOL = set_a(obj,value,varargin{:});
                case 'eff'
                    BOOL = set_eff(obj,value,varargin{:});
                case 'pf'
                    BOOL = set_pf(obj,value,varargin{:});
                case 'Span'
                    BOOL = set_Span(obj,value,varargin{:});
                case 'Lend'
                    BOOL = set_Lend(obj,value,varargin{:});
                case 'thetar'
                    BOOL = set_thetar(obj,value,varargin{:});
                case 'thetag'
                    BOOL = set_thetag(obj,value,varargin{:});
                case 'thetas'
                    BOOL = set_thetas(obj,value,varargin{:});
                case 'Bg'
                    BOOL = set_Bg(obj,value,varargin{:});
                case 'Bt'
                    BOOL = set_Bt(obj,value,varargin{:});
                case 'Bry'
                    BOOL = set_Bry(obj,value,varargin{:});
                case 'Bsy'
                    BOOL = set_Bsy(obj,value,varargin{:});
                case 'winT'
                    BOOL = set_winT(obj,value,varargin{:});
                case 'alphaT'
                    BOOL = set_alphaT(obj,value,varargin{:});
                case 'kR'
                    BOOL = set_kR(obj,value,varargin{:});
            end
        end
        function BOOL = set_Pout(obj,value,varargin)
            if value<0
                PrintError('Output power can not be a negative number',varargin{:});
                BOOL = true;
            else
                obj.Pout = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_rpm(obj,value,varargin)
            if value<0
                PrintError('Rated Speed can not be a negative number',varargin{:});
                BOOL = true;
            else
                obj.rpm = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Vrms(obj,value,varargin)
            if value<0
                PrintError('RMS of phase voltage can not be a negative number',varargin{:});
                BOOL = true;
            else
                obj.Vrms = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_kR(obj,value,varargin)
            if value<1
                PrintError('kR must be grater than 1',varargin{:});
                BOOL = true;
            else
                obj.kR = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Ns(obj, value,varargin)
            if value<0 || rem(value,1)
                PrintError('Number of stator slot can not be a negative number',varargin{:});
                BOOL = true;
            else
                obj.Ns = value;
                obj.isSetGeom = false;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Nm(obj, value,varargin)
            if value<0
                PrintError('Number of magnets can not be a negative number',varargin{:});
                BOOL = true;
            elseif rem(value,2)
                PrintError('Number of magnets must be an even integer',varargin{:});
                BOOL = true;
            else
                obj.Nm = value;
                obj.isSetGeom = false;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_g(obj,value,varargin)
            if value<0
                PrintError('Air gap length must be a positive number',varargin{:});
                BOOL = true;
            else
                obj.g = value;
                obj.isSetGeom = false;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Rso(obj,value,varargin)
            if value<0
                PrintError('Outer stator radius must be a positive number',varargin{:});
                BOOL = true;
            else
                obj.Rso = value;
                obj.isSetGeom = false;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_hm(obj,value,varargin)
            if value<0
                PrintError('Magnet height length must be a positive number',varargin{:});
                BOOL = true;
            else
                obj.hm = value;
                obj.isSetGeom = false;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_winT(obj,value,varargin)
            if value<0
                PrintError('Winding Temperature must be a positive value',varargin{:});
                BOOL = true;
            else
                obj.winT = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_alphaT(obj,value,varargin)
            if value<0
                PrintError('alphaT must be a positive value',varargin{:});
                BOOL = true;
            else
                obj.alphaT = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Lst(obj,value,varargin)
            if value<0
                PrintError('Stack length can not be a negative number',varargin{:});
                BOOL = true;
            else
                obj.Lst = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Ncoil(obj, value,varargin)
            if value<0
                PrintError('Number of coil turns can not be a negative number',varargin{:});
                BOOL = true;
            elseif rem(value,1)
                PrintError('Number of coil turns must be a positive integer',varargin{:});
                BOOL = true;
            else
                obj.Ncoil = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_betam(obj, value,varargin)
            if value<0.5 || value > 0.95
                PrintError('betam must be between 0.5 - 0.95',varargin{:});
                BOOL = true;
            else
                obj.betam = value;
                obj.isSetGeom = false;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_m(obj, value,varargin)
            if value<0 || rem(value,1)
                PrintError('Number of phases must be a positive integer',varargin{:});
                BOOL = true;
            else
                obj.m = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Kst(obj, value,varargin)
            if value<0.9 || value > 1
                PrintError('Staking factor must be between 0.9 - 1',varargin{:});
                BOOL = true;
            else
                obj.Kst = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_D(obj, value,varargin)
            if value<0.45 || value > 0.7
                PrintError('betam must be between 0.45 - 0.7',varargin{:});
                BOOL = true;
            else
                obj.D = value;
                obj.isSetGeom = false;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Kf(obj, value,varargin)
            if value<0.2 || value > 0.7
                PrintError('Winding fill factor must be between 0.2 - 0.7',varargin{:});
                BOOL = true;
            else
                obj.Kf = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Bg(obj, value,varargin)
            if value<0.5 || value > 1.5
                PrintError('Bg must be between 0.5 - 1.5',varargin{:});
                BOOL = true;
            else
                obj.Bg = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Bt(obj, value,varargin)
            if value<1 || value > 2.5
                PrintError('Bt must be between 1 - 2.5',varargin{:});
                BOOL = true;
            else
                obj.Bt = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Bsy(obj, value,varargin)
            if value<1 || value > 2
                PrintError('Bsy must be between 1 - 2',varargin{:});
                BOOL = true;
            else
                obj.Bsy = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Bry(obj, value,varargin)
            if value<1 || value > 2
                PrintError('Bry must be between 1 - 2',varargin{:});
                BOOL = true;
            else
                obj.Bry = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_eff(obj, value,varargin)
            if value<0.4 || value > 1
                PrintError('Efficiency must be between 0.4 - 1',varargin{:});
                BOOL = true;
            else
                obj.eff = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_pf(obj, value,varargin)
            if value<0.5 || value > 1
                PrintError('Power Factor must be between 0.5 - 1',varargin{:});
                BOOL = true;
            else
                obj.pf = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_wry(obj, value,varargin)
            if value<0
                PrintError('Number of phases must be a positive number',varargin{:});
                BOOL = true;
            else
                obj.wry = value;
                obj.isSetGeom = false;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_wsy(obj, value,varargin)
            if value<0
                PrintError('Number of phases must be a positive number',varargin{:});
                BOOL = true;
            else
                obj.wsy = value;
                obj.isSetGeom = false;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_wtb(obj, value,varargin)
            if value<0
                PrintError('Number of phases must be a positive number',varargin{:});
                BOOL = true;
            else
                obj.wtb = value;
                obj.isSetGeom = false;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_bss1(obj, value,varargin)
            if value<0
                PrintError('Number of phases must be a positive number',varargin{:});
                BOOL = true;
            else
                obj.bss1 = value;
                obj.isSetGeom = false;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_hss1(obj, value,varargin)
            if value<0
                PrintError('Number of phases must be a positive number',varargin{:});
                BOOL = true;
            else
                obj.hss1 = value;
                obj.isSetGeom = false;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_a(obj, value,varargin)
            if value<0 || rem(value,1)
                PrintError('Number of parallel path must be a positive number',varargin{:});
                BOOL = true;
            else
                obj.a = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Span(obj, value,varargin)
            if value<0 || rem(value,1)
                PrintError('Coil span must be a positive integer',varargin{:});
                BOOL = true;
            else
                obj.Span = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_Lend(obj, value,varargin)
            if value<0
                PrintError('Lend must be positive number',varargin{:});
                BOOL = true;
            else
                obj.Lend = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_thetar(obj,value,varargin)
            if value<0
                PrintError('thetar must be a positive number',varargin{:});
                BOOL = true;
            else
                obj.thetar = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_thetas(obj,value,varargin)
            if value<0
                PrintError('thetas must be a positive number',varargin{:});
                BOOL = true;
            else
                obj.thetas = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = set_thetag(obj,value,varargin)
            if value<0
                PrintError('thetag must be a positive number',varargin{:});
                BOOL = true;
            else
                obj.thetag = value;
                BOOL = false;
                PrintGoodMessage(varargin{:});
            end
        end
        function BOOL = checkIterdependencies(obj,varargin)
            if rem(obj.Ns/obj.m,1)
                PrintError('the ratio of Ns/m must be a positive integer, number of coils in each phase must be equal',varargin{:});
            elseif rem(obj.Ns/obj.m/obj.a,1)
                PrintError('Number of coils in each parallel path must be equal, you must change Ns or m or a',varargin{:});
            else
                BOOL = false;
                return;
            end
            BOOL = true;
        end
        %% Evaluation of dependent variables
        function y = get.wm(obj)
            y = obj.rpm*pi/30;
        end
        function y = get.we(obj)
            y = obj.wm*obj.Nm/2;
        end
        function y = get.T(obj)
            y = obj.Pout*1e3/obj.wm;
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
            y = 2*asin(obj.bss1/2/(obj.Rro+obj.g));
        end
        function y = get.rss1(obj)
            ttmp = asin(obj.bss1/2/(obj.Rro+obj.g));
            Atmp = [1,1;1,sin(obj.Ts/2)];
            btmp = [obj.Rso-obj.wsy-obj.rss2-(obj.Rro+obj.g)*cos(ttmp)-obj.hss1-...
                tan(obj.Ts/2)*obj.bss1/2;
                (obj.Rso-obj.wsy-obj.rss2)*sin(obj.Ts/2)-obj.wtb/2-obj.bss1/2/cos(obj.Ts/2)];
            xtmp = Atmp\btmp;
            y = xtmp(1);
        end
        function y = get.rss2(obj)
            y = ((obj.Rso-obj.wsy)*sin(obj.Ts/2)-obj.wtb/2)/(1+sin(obj.Ts/2));
        end
        function y = get.f(obj)
            y = obj.rpm*obj.Nm/120;
        end
        function y = get.Rdc_coil(obj)
            Acopper = (obj.Aslot/2)*obj.Kf;
            Rslot = obj.Ncoil^2*obj.Lst/(Acopper*obj.sigma*1e-3);
            Rave = 0.5*(obj.Rro+obj.g+obj.Rso);
            L2 = Rave*obj.Span*obj.Ts;
            Lew = L2 + 2*obj.Lend;
            Rew = obj.Ncoil^2*Lew/(Acopper*obj.sigma*1e-3);
            y = 2*(Rslot+Rew);
        end
        function y = get.Rdc_coil_T(obj)
            y = obj.Rdc_coil*(1+obj.alphaT*(obj.winT-25));
        end
        function y = get.Rac_coil_T(obj)
            y = obj.Rdc_coil_T*obj.kR;
        end
        function y = get.Rph(obj)
            y = obj.Ns*obj.Rac_coil_T/obj.a^2/obj.m;
        end
        function y = get.Rro(obj)
            y = obj.D*obj.Rso;
        end
        function y = get.CopperLoss(obj)
            y = obj.m*obj.Rph*obj.Irms^2;
        end
        function y = get.Nph(obj)
            y = obj.Ns*obj.Ncoil/(obj.a*obj.m);
        end
        function y = get.Irms(obj)
            Pin = obj.Pout*1e3/obj.eff;
            y = Pin/obj.m/obj.Vrms/obj.pf;
        end
        function y = get.phi_tot(obj)
            y = 2*pi*obj.Rro*obj.Bg*obj.Lst*1e-6;
        end
        function y = get.phi_p(obj)
            y = obj.phi_tot/obj.Nm;
        end
        function update_wtb_For_Bt(obj)
            obj.wtb = 2*pi*obj.Rro*obj.Bg/(obj.Ns*obj.Kst*obj.Bt);
            obj.wtb = round(obj.wtb,1);
        end
        function update_wsy_For_Bsy(obj)
            obj.wsy = pi*obj.Rro*obj.Bg/(obj.Nm*obj.Kst*obj.Bsy);
            obj.wsy = round(obj.wsy,1);
        end
        function update_wry_For_Bry(obj)
            obj.wry = pi*obj.Rro*obj.Bg/(obj.Nm*obj.Kst*obj.Bry);
            obj.wry = round(obj.wry,1);
        end
        function update_Ncoil_For_BEMF(obj)
            obj.Ncoil = obj.Vrms/(4.44*obj.f*obj.phi_p*0.955);
            obj.Ncoil = obj.a*obj.m*obj.Ncoil/obj.Ns;
            obj.Ncoil = ceil(obj.Ncoil);
        end
        function update_geom(obj)
            if ~obj.isSetGeom
                obj.writeParFile;
                obj.geom.clearAllmzs;
                geom_pmsm0;
                obj.geom.read_g2d_bin('geom.g2d','MM');
                obj.isSetGeom = true;
            end   
        end
        function writeParFile(obj)
            writeParFile('Rso',obj.Rso,'hm',obj.hm,'Rro',obj.Rro,...
                'rsh',obj.rsh,'Ts',obj.Ts,'Tr',obj.Tr,'thetag',obj.thetag,...
                'thetar',obj.thetar,'thetas',obj.thetas,...
                'bss1',obj.bss1,'hss1',obj.hss1,'rss1',obj.rss1,...
                'rss2',obj.rss2,'Tso1',obj.Tso1,...
                'wry',obj.wry,'wsy',obj.wsy,'g',obj.g,'betam',obj.betam);
        end
        function y = get.Aslot(obj)
            if ~obj.isSetGeom
                obj.update_geom;
            end
            y = 2*obj.geom.mzs.c11.area;
        end
        function obj = setWin(obj,Directory)
            File = fopen(Directory,'r');
            if File<0
                error('program can not find file\n');
            end
            tmp = fscanf(File,'%d');
            fclose(File);
            % check for consistensy with the number of slots
            if obj.Ns ~= tmp(1)
                error(['Windings data does not set properly',...
                    ' ... Number of coils must be equal to the number of slots\n']);
            end
            % check for consistensy with the number of phases
            if tmp(2) ~= obj.m
                error('Phase number is wrong ...\n');
            end
            % check for double layer winding
            if tmp(4) ~= 2
                error('Number of layers should be 2 ...\n');
            end
            obj.win.Span = tmp(3);
            % getting phases characteristics
            tmp = reshape(tmp(5:end),2,[]);
            tmp = tmp';
            tmp = reshape(tmp,[],2*obj.m);
            for i = 1:obj.m
                phase = tmp(:,[i,i+obj.m]);
                if any(bitor(phase(:,1)>obj.Ns,phase(:,1)<1))
                    error(['Windings data does not set properly',...
                        ' ... Some coil index is grater than Ns\n'])
                end
                obj.win.(['Phase',num2str(i)]) = phase;
            end
            % setting coils
            obj.win.Coils = (1:obj.Ns)';
            obj.win.Coils = [obj.win.Coils,circshift(obj.win.Coils,-obj.win.Span)];
            % set flag
            obj.isSetWin = true;
        end
        function readFile(obj,Dir)
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
                if isfield(obj,str{1})
                    obj.(str{1}) = str2double(str{2});
                end
            end
        end
        function writeFile(obj,Dir)
            if nargin<2
                Dir = cd;
            end
            File = fopen([Dir,'\data.par'],'w');
            if ~File
                error('Wrong file directory');
            end
            varNames = fieldnames(obj);
            for i = 1:numel(varNames)
                if ~isstruct(obj.(varNames{i}))
                    fprintf(File,'%s\t%20.20f\n',varNames{i},obj.(varNames{i}));
                end
            end
            fclose(File);
            copyfile data.pmsm G2DKERNEL\data.par
        end
    end
end
function PrintError(Message,varargin)
if nargin == 1
    error(Message);
elseif nargin == 2
    varargin{1}.String = '';
    varargin{1}.ForegroundColor = 'r';
    varargin{1}.String{1} = 'Error:';
    varargin{1}.String{2} = [Message,' ...'];
else
    error('INTERNAL ERROR');
end
end
function PrintGoodMessage(varargin)
if nargin == 1
    varargin{1}.String = '';
    varargin{1}.ForegroundColor = 'b';
    varargin{1}.String{1} = 'Done!';
    varargin{1}.String{2} = 'Parameter has been set correctly.';
elseif nargin>1
    error('INTERNAL ERROR');
end
end