classdef emdlab_mlib_electricalSteel < emdlab_phy_material & handle
    
    properties
        
        gradeName (1,:) char;
        
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
        
        function obj = emdlab_mlib_electricalSteel()
            
            obj.MagneticPermeability.isIsotropic = true;
            obj.MagneticPermeability.isLinear = false;
            
        end
        
        function plotHBCurve(obj, aHandle)
            
            if nargin < 2
                aHandle = gca;
            end
            
            x = obj.hb_curve(:,1)/1000;
            plot(aHandle, x, obj.hb_curve(:,2), 'LineWidth',1.5);
            box(aHandle, 'on');
            xlabel(aHandle, 'H [kA/m]');
            ylabel(aHandle, 'B [tesla]');
            title(aHandle, obj.gradeName);
            grid(aHandle, 'on');
            aHandle.XAxis.FontSize = 12;
            aHandle.YAxis.FontSize = 12;
            aHandle.Title.FontSize = 14;
            aHandle.XLim(1) = -max(x)/10;
            
        end
        
        function plotHBCurveLog(obj, aHandle)
            
            if nargin < 2
                aHandle = gca;
            end
            
            plot(aHandle, log(obj.hb_curve(:,1)), obj.hb_curve(:,2));
            box(aHandle, 'on');
            xlabel(aHandle, 'log(H) [A/m]');
            ylabel(aHandle, 'B [Tesla]');
            title(aHandle, obj.gradeName);
            grid(aHandle, 'on');
            
        end
        
        function evalHBCurveRelatedQuantities(obj)
            
            h = obj.hb_curve(2:end,1);
            b = obj.hb_curve(2:end,2);

            % smooth & extend HB curve
            [h,b] = emdlab_flib_smooth_extend_hbcurve_exp(h, b, 100*h(end));
            %[h,b] = emdlab_flib_smooth_extend_hbcurve_polyc(h, b, 10*h(end));
            %[h,b] = emdlab_flib_smooth_extend_hbcurve_arctan(h, b, 10*h(end));

            % evaluation of nu
            v = (h ./ b);
%             [b, v] = extend_exponential(b, v, 1/4/pi/1e-7);
            %[b, v] = emdlab_flib_smooth_extend_bvcurve_arctan(b,v,1/4/pi/1e-7);
            %[b, v] = emdlab_flib_smooth_extend_bvcurve_exp(b,v,1/4/pi/1e-7);
            h = v .* b;
            
            % evaluation of BH, vB and dvdB
            h = [-flipud(h(2:end)); 0; h];
            b = [-flipud(b(2:end)); 0; b];
            v = [flipud(v(2:end)); v(1); v];
            
            obj.BH = pchip(b, h);
            obj.dBdH = obj.BH;
            obj.dBdH.coefs = obj.dBdH.coefs * diag(3:-1:1, 1);
            
            obj.vB = pchip(b, v);
            obj.dvdB = obj.vB;
            obj.dvdB.coefs = obj.dvdB.coefs * diag(3:-1:1, 1);
            
%             Ntmp = 10000;
%             obj.be = linspace(0, 10, Ntmp)';
%             beStep = obj.be(2);
%             beMid = (obj.be(1:end-1) + obj.be(2:end))/2;
%             vMid = ppval(obj.vB, beMid);
%             obj.we = cumsum([0; vMid.*beMid*beStep]);

            Ntmp = 1000;
            obj.be = linspace(0, 10, Ntmp)';
            nuPoints = ppval(obj.BH, obj.be);
            obj.we = zeros(1,Ntmp);
            for i = 2:Ntmp                
                obj.we(i) = trapz(obj.be(1:i), nuPoints(1:i));
            end
            
        end
        
        function writeMaxwellScript(obj, e2m)
            
            e2m.writeLine('Set oDefinitionManager = oProject.GetDefinitionManager()');
            e2m.writeLine(sprintf('oDefinitionManager.AddMaterial Array("NAME:%s", "CoordinateSystemType:=",  _', obj.gradeName));
            e2m.writeLine('"Cartesian", Array("NAME:AttachedData"), Array("NAME:ModifierData"), Array("NAME:permeability", "property_type:=",  _');
            e2m.writeLine('"nonlinear", "BType:=", "normal", "HUnit:=", "A_per_meter", "BUnit:=", "tesla", Array("NAME:BHCoordinates", _');
            
            for i = 1:length(obj.hb_curve)-1
                
                e2m.writeLine(sprintf('Array("NAME:Coordinate", "X:=", %f, "Y:=", %f), _', obj.hb_curve(i,1), obj.hb_curve(i,2)));
                
            end
            e2m.writeLine(sprintf('Array("NAME:Coordinate", "X:=", %f, "Y:=", %f))), _', obj.hb_curve(end,1), obj.hb_curve(end,2)));
            
            e2m.writeLine('Array("NAME:magnetic_coercivity", "property_type:=",  _');
            e2m.writeLine('"VectorProperty", "Magnitude:=", "0A_per_meter", "DirComp1:=", "1", "DirComp2:=",  _');
            e2m.writeLine('"0", "DirComp3:=", "0"), Array("NAME:stacking_type", "property_type:=",  _');
            e2m.writeLine('"ChoiceProperty", "Choice:=", "Lamination"), "stacking_factor:=", "0.95", Array("NAME:stacking_direction", "property_type:=",  _');
            e2m.writeLine('"ChoiceProperty", "Choice:=", "V(3)"))');   
            
        end
        
    end
    
    
end