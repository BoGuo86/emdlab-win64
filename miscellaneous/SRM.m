classdef SRM
    properties
        
        Pout
        rpm
             
        Rso = 50;
        D = 0.55;
        
        betas
        betar
        
        Nr = 8;
        Ns = 12;
        
        wsyPwsp = 0.6;
        wryPwrp = 0.6;
        
        
        
    end
    properties (Dependent = true)
        
        wm
        T
        
        Rro
        wsp
        wrp
        wry
        wsy
        Nsplit
        
        
    end
    methods
        function y = get.wm(obj)
            y = obj.rpm*pi/30;
        end
         function y = get.T(obj)
            y = obj.Pout/obj.wm;
         end
        
         function y = get.Nsplit(obj)
             y = gcd(obj.Ns,obj.Nr);
         end
    end
end