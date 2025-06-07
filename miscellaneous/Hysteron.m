classdef Hysteron < handle
    
    properties
        
        output (1,1) double;
        yLow (1,1) double = 0;
        yHigh (1,1) double = 1;
        xLow (1,1) double = -0.5;
        xHigh (1,1) double = 0.5;
        
    end
    
    methods
        
        function obj = Hysteron(xLow, xHigh, yLow, yHigh)
            
            obj.xLow = xLow;
            obj.xHigh = xHigh;
            obj.yLow = yLow;
            obj.yHigh = yHigh;
            obj.output = yLow;
            
        end
        
        function updateOutput(obj, input)
            
            if input >= obj.xHigh
                
                obj.output = obj.yHigh;
                
            elseif input <= obj.xLow
                
                obj.output = obj.yLow;
                
            end
            
        end
        
        function y = getOutput(obj)

            y = obj.output;
            
        end
        
        function setOutputHigh(obj)
            
            obj.output = obj.yHigh;
            
        end
        
        function setOutputLow(obj)
            
            obj.output = obj.setOutputLow;
            
        end
        
    end
    
end
