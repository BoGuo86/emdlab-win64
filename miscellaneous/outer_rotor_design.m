classdef outer_rotor_design < handle
    properties
            
        g
        hm
        wry
        wsy
        wtb
        
        bss1
        hss1
        Tsso1
        Tsso2
        Tsso3
        
        Rri
        Rro
        
        Rsi
        Rso
        
        D
        
        Ns
        Nm
        
        Ts
        Tr
        
        rpm
        rps
        wm
        we
        
        Lst
        
        Nsplit
        m
        
    end
    methods
        function obj = outer_rotor_design()
        end
        function set.Rro(obj,value)
            if value<0
                error()
            else
                obj.Rro = value;
            end
        end
        function y = get.Rso(obj)
            y = obj.D*obj.Rro;
        end
    end
end