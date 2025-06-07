classdef emdlab_winding_3ph2layer < handle
    
    properties
        
        % winding data
        Ns (1,1) double;
        p (1,1) double;
        span (1,1) double;
        
        % electrical angle & mechanical angle of coil
        eaCoil (:,1) double;
        maCoil (:,1) double;
        
        % coil span angle
        eaSpan
        maSpan
        
        % harmonic content of generated mmf by winding
        an
        bn
        
        coils
        phases
        
        % mechanical angle of magnetic axis of phase a
        mamaa
        
        slot2PhaseName
        slot2CoilSign
        
    end
    
    methods
        
        function obj = emdlab_winding_3ph2layer(varargin)
            
            if numel(varargin) == 0
                
                obj.Ns = 36;
                obj.p = 4;
                obj.span = 5;
                
            elseif numel(varargin) == 2
                
                obj.Ns = varargin{1};
                obj.p = varargin{2};
                obj.span = max(floor(obj.Ns/obj.p-1), 1);
                
            elseif numel(varargin) == 3
                
                obj.Ns = varargin{1};
                obj.p = varargin{2};
                obj.span = varargin{3};
                
            else
                
                MException('', 'Incorrect number of inputs.');
                
            end
            
%             obj.updateWinding;       
            
        end
        
        function setWinding(obj, Ns_, p_, span_)
            
            obj.Ns = Ns_;
            obj.p = p_;
            obj.span = span_;
            obj.updateWinding;
            
        end
        
        function updateWinding(obj)

            obj.coils = (1:obj.Ns)';
            obj.coils = [obj.coils, circshift(obj.coils, -obj.span)];
            
            slot_angle = 2*pi/obj.Ns;
            obj.maSpan = obj.span*slot_angle;
            obj.maCoil = (0:obj.Ns-1) * slot_angle + obj.maSpan/2;
            
            alpha = (0:obj.Ns-1) * (obj.p/2) * (2*pi/obj.Ns);
            alpha = rem(alpha, 2*pi-1e-6);
            index = alpha>=pi;
            alpha(index) = alpha(index)-pi;
            pts = exp(complex(0,1)*alpha);
            
            tmp = 1:obj.Ns;
            tmp(index) = -tmp(index);
            [~, index2] = sort(abs(pts-pts(1)));
            tmp = tmp(index2);
            obj.phases = [tmp(1:obj.Ns/3);tmp(2*obj.Ns/3+1:obj.Ns);-tmp(obj.Ns/3+1:2*obj.Ns/3)]';
            
            [~,index] = sort(abs(obj.phases));
            obj.phases(:,1) = obj.phases(index(:,1),1);
            obj.phases(:,2) = obj.phases(index(:,2),2);
            obj.phases(:,3) = obj.phases(index(:,3),3);
       
            phaseNames = {'phaseA', 'phaseB', 'phaseC'};
            
            obj.slot2PhaseName = cell(obj.Ns,2);
            obj.slot2CoilSign = zeros(obj.Ns,2);
            
            for j = 1:3
                
                for i = 1:obj.Ns/3
                    
                    cIndex = obj.coils(abs(obj.phases(i,j)),:);
                    obj.slot2PhaseName{cIndex(1),1} = phaseNames{j};
                    obj.slot2PhaseName{cIndex(2),2} = phaseNames{j};
                    
                    if obj.phases(i,j) > 0
                        
                        obj.slot2CoilSign(cIndex(1),1) = 1;
                        obj.slot2CoilSign(cIndex(2),2) = -1;
                        
                    else
                        
                        obj.slot2CoilSign(cIndex(1),1) = -1;
                        obj.slot2CoilSign(cIndex(2),2) = 1;
                        
                    end
                    
                end
                
            end
            
        end
        
        function plotMMFDistribution(obj)
            
            figure;
            subplot(211);
            hold on;
            
            t = linspace(0,2*pi,1000);
            obj.eval_an_bn;
            
            v = 1:300;
            
            Ncpph = obj.Ns/obj.p/3;
            
            y = obj.an(1,:)*cos(v'*t)+obj.bn(1,:)*sin(v'*t);
            yt = y;
            plot(t*180/pi, y/Ncpph, 'LineWidth', 1.5, 'color', 'b');
            
            y = obj.an(2,:)*cos(v'*t)+obj.bn(2,:)*sin(v'*t);
            yt = yt -0.5*y;
            
            plot(t*180/pi, y/Ncpph, '--', 'color', 'r');
            
            y = obj.an(3,:)*cos(v'*t)+obj.bn(3,:)*sin(v'*t);
            yt = yt -0.5*y;
            plot(t*180/pi, y/Ncpph, '--', 'color', 'k');

            set(gca,'xlim', [0,360]);
            
            
            h = sqrt(obj.an(1,:).^2+obj.bn(1,:).^2);
            [~,index] = max(h);
            y = obj.an(1,index)*cos(index*t)+obj.bn(1,index)*sin(index*t);
            plot(t*180/pi, y/Ncpph, '--', 'color', 'b');
            
            [~, tmp] = sort(abs(y));
            obj.mamaa = sort(t(tmp(end-obj.p:end)));
            
            h1 = h(index);
            h(index) = 0;
            thd = sqrt(sum(h.^2))*100/h1;
            
            str_title = sprintf('Kw = %.3f, THD = %.2f%%', max(y/Ncpph)/(4/pi), thd);
            
            title(str_title);
            legend({'Phase A', 'Phase B', 'Phase C', 'First Order'});
            
            box('on');
            xlabel('Mechanical Angle [Deg]');
            ylabel('MMF (per unit) [A]');
            
subplot(212);
            plot(t*180/pi, yt/1.5/Ncpph, 'LineWidth', 1.5, 'color', 'b');

            
            tmobj.pan = obj.an(1,:)-0.5*obj.an(2,:)-0.5*obj.an(3,:);
            tmobj.pbn = obj.bn(1,:)-0.5*obj.bn(2,:)-0.5*obj.bn(3,:);
            h = sqrt(tmobj.pan.^2+tmobj.pbn.^2);
            [~,index] = max(h);
            y = tmobj.pan(index)*cos(index*t)+tmobj.pbn(index)*sin(index*t);
            plot(t*180/pi, y/Ncpph/1.5, '--', 'color', 'b');
            
            h1 = h(index);
            h(index) = 0;
            thd = sqrt(sum(h.^2))*100/h1;
            
            str_title = sprintf('Kw = %.3f, THD = %.2f%%', max(y/Ncpph/1.5)/(4/pi), thd);
            title(str_title);
            
            box('on');
            xlabel('Mechanical Angle [Deg]');
            ylabel('MMF (per unit) [A]');
            set(gca, 'xlim', [0,360]);
            
        end
        
        function eval_an_bn(obj)
            
            v = 1:300;
            obj.an = zeros(3,300);
            obj.bn = zeros(3,300);
            tmp = (2/pi./v).*sin(v*obj.maSpan/2);
            
            
            for i = 1:3
                
                for j = 1:obj.Ns/3
                    
                    cIndex = obj.phases(j,i);
                    obj.an(i,:) = obj.an(i,:) + sign(cIndex) * cos(v*obj.maCoil(abs(cIndex)));
                    obj.bn(i,:) = obj.bn(i,:) + sign(cIndex) * sin(v*obj.maCoil(abs(cIndex)));
                    
                end
                
                obj.an(i,:) = tmp.*obj.an(i,:);
                obj.bn(i,:) = tmp.*obj.bn(i,:);
                
            end
            
            
        end
        
    end
    
end