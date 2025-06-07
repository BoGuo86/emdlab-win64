classdef ID_IM1PHASE0
    
    properties
        POUT = 100;   
        RPM = 1400;
        P = 4;
        Ns = 16;
        
        RSO = 100;
        D = 0.6;  
    end
    
    properties(Dependent = true)
        RPS
        TOUT
        WM
        WE
        RRO
    end
    
    methods
        
        function obj = ID_IM1PHASE0()
        end
                      
        function y = get.WM(obj)
            y = obj.RPM * pi/30;
        end
        
        function y = get.WE(obj)
            y = obj.WM * obj.P/2;
        end
        
        
    end
    
end