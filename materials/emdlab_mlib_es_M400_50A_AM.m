classdef emdlab_mlib_es_M400_50A_AM < handle & emdlab_mlib_electricalSteel
    
    
    methods
        
        function obj = emdlab_mlib_es_M400_50A_AM()
            
            obj.gradeName = 'M400-50A-AM';
            
            obj.hb_curve = [0	0
5	0.0081679
10	0.0174667
20	0.044735
25	0.06283
32.6	0.1
43.5	0.2
50.8	0.3
57.2	0.4
63.4	0.5
69.9	0.6
77.3	0.7
86	0.8
97.2	0.9
113.2	1
137.8	1.1
180.2	1.2
269.5	1.3
516.8	1.4
1307	1.5
3180	1.6
6361	1.7
10890	1.79526
20000	1.8935
40000	2.00551
80000	2.11752
120000	2.18304
140000	2.20795];
            
            obj.evalHBCurveRelatedQuantities;
            
        end
        
    end
    
end