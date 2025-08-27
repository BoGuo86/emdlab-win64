classdef emdlab_mlib_permanentMagnet < emdlab_phy_material & handle
    
    properties
        
        
        
        hb_curve (:,2) double;
        BH
        dBdH
        vB
        vH
        dvdB
        vB2
        dvdB2
        be
        we
        wc
        
    end
    
    methods
        
        function obj = emdlab_mlib_permanentMagnet()
            
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isLinear = true;
            
        end
        
        function plotHBCurve(obj, aHandle)
            
            if nargin < 2
                aHandle = gca;
            end
            
            plot(aHandle, obj.hb_curve(:,1), obj.hb_curve(:,2));
            box(aHandle, 'on');
            xlabel(aHandle, 'H [A/m]');
            ylabel(aHandle, 'B [Tesla]');
            grid(aHandle, 'on');
            
        end
        
        function plotHBCurveLog(obj, aHandle)
            
            if nargin < 2
                aHandle = gca;
            end
            
            plot(aHandle, log(obj.hb_curve(:,1)), obj.hb_curve(:,2));
            box(aHandle, 'on');
            xlabel(aHandle, 'log(H) [A/m]');
            ylabel(aHandle, 'B [Tesla]');
            grid(aHandle, 'on');
            
        end
        
        function evalHBCurveRelatedQuantities(obj)
            
            % getting and smoothing HB curve
            [h, b] = smoothHBCurve(obj.hb_curve(2:end,:), 135);
            
            % evaluation of nu
            v = (h ./ b);
            [b, v] = extend_exponential(b, v, 1/4/pi/1e-7);
            h = v .* b;
            
            % evaluation of BH, vB and dvdB
            h = [-flipud(h(2:end)); 0; h];
            b = [-flipud(b(2:end)); 0; b];
            v = [flipud(v(2:end)); v(1); v];
            
            obj.BH = spline(b, h);
            obj.dBdH = obj.BH;
            obj.dBdH.coefs = obj.dBdH.coefs * diag(3:-1:1, 1);
            
            obj.vB = spline(b, v);
            obj.dvdB = obj.vB;
            obj.dvdB.coefs = obj.dvdB.coefs * diag(3:-1:1, 1);
            
            Ntmp = 1000;
            obj.be = linspace(0, 5, Ntmp)';
            beStep = obj.be(2);
            beMid = (obj.be(1:end-1) + obj.be(2:end))/2;
            vMid = ppval(obj.vB, beMid);
            obj.we = cumsum([0; vMid.*beMid*beStep]);
            
        end
        
    end
    
    
end