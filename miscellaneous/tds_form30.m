
classdef tds_form30 < handle
   
    properties
        DS (1,1) double {mustBePositive} = 1090;
        HF (1,1) double {mustBePositive} = 2248;
        ES (1,1) double {mustBePositive} = 2019;
        BS (1,1) double {mustBePositive} = 1000;
        SS (1,1) double {mustBePositive} = 923;
        BJ (1,1) double {mustBePositive} = 1000;
    end
    
    properties
        p (1,1) matlab.graphics.primitive.Patch
    end
    
    methods
        function obj = tds_form30(f)
            obj.p.Parent = (f);
        end
    end
    
end