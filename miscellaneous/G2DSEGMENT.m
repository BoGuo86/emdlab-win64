% line properties class
classdef G2DSEGMENT
    properties (Constant)
        % geometry epsilon for length
        gleps = 1e-6 
        % geometry epsilon for angle
        gaeps = 2*atan(1e-6/2) 
    end 
    properties
        p1
        p2
        maxLength
        Nnodes
        alpha
        L1
        L2
    end
    methods  
        function obj = G2DSEGMENT(p1,p2,varargin)
            obj.p1 = p1;
            obj.p2 = p2;
            
            p = inputParser();
            addParameter(p,'maxLength',inf);
            addParameter(p,'Nnodes',2);
            addParameter(p,'alpha',inf);
            addParameter(p,'L1',inf);
            addParameter(p,'L2',inf);
            p.CaseSensitive = 0;
            parse(p,varargin{:});
            
            obj.maxLength = p.Results.maxLength;
            obj.Nnodes = p.Results.Nnodes;
            obj.alpha = p.Results.alpha;
            obj.L1 = p.Results.L1;
            obj.L2 = p.Results.L2;

        end
        function y = getip(obj)
            L = obj.getLength;
            u = obj.getu;
            if obj.Nnodes>2
                % check for alpha
                if obj.alpha ~= inf
                    if obj.alpha ~= 1
                        d = L * (1-obj.alpha)/(1-obj.alpha^(obj.Nnodes-1));
                        r = obj.alpha;
                    end
                elseif obj.L1 ~= inf
                    if obj.L2 ~= inf
                        d = obj.L1;
                        r = nthroot(obj.L2/obj.L1,obj.Nnodes-2);
                    end
                else
                    d = L/(obj.Nnodes-1);
                        r = 1;
                    
                end
                
                
                
                y = zeros(obj.Nnodes,2);
                y(1,:) = obj.p1;
                for i = 2:obj.Nnodes
                    y(i,:) = y(i-1,:) + u*d*r^(i-2);
                end
                y = y(2:end-1,:);
                
            elseif obj.L1 ~= inf
                if obj.L2 ~= inf
                    if abs(obj.L1-obj.L2)<obj.gleps
                        n = L/obj.L1;
                        r = 1;
                    else
                        r = (L-obj.L1)/(L-obj.L2);
                        n = round(log10(obj.L2/obj.L1)/log10(r))+1;
%                         r = nthroot(obj.L2/obj.L1,n-1);
                    end 
                else
                    n = L/obj.L1;
                    r = 1; 
                end
                y = zeros(n+1,2);
                y(1,:) = obj.p1;
                for i = 2:n+1
                    y(i,:) = y(i-1,:) + u*obj.L1*r^(i-2);
                end
                y = y(2:end-1,:);
            else
                if obj.maxLength == inf
                    y = zeros([],2);
                else
                    L = floor(obj.getLength/obj.gleps)*obj.gleps;
                    mL = floor(obj.maxLength/obj.gleps)*obj.gleps;
                    n = ceil(L/mL)+1;
                    y = [linspace(obj.p1(1),obj.p2(1),n);...
                        linspace(obj.p1(2),obj.p2(2),n)]';
                    y = y(2:end-1,:);
                end
            end
        end
        function [e,v] = getev(obj)
            v = [obj.p1;obj.getip;obj.p2];
            Nv = size(v,1);
            e = [1:Nv-1;2:Nv]';
        end
        function y = getLength(obj)
            y = norm(obj.p2-obj.p1);
        end
        function y = getu(obj)
            y = (obj.p2-obj.p1)/obj.getLength;
        end
        function y = getmirror(obj,varargin)
            y = obj;
            y.p1 = pmirror(y.p1,varargin{:});
            y.p2 = pmirror(y.p2,varargin{:});
        end
        function y = getrotate(obj,varargin)
            y = obj;
            y.p1 = ext_protate2(y.p1,varargin{:});
            y.p2 = ext_protate2(y.p2,varargin{:});
        end
        function y = getshift(obj,varargin)
            y = obj;
            y.p1 = pshift(y.p1,varargin{:});
            y.p2 = pshift(y.p2,varargin{:});
        end
        function y = getindent(obj,d)
            y = obj;
            u = y.getu;
            u = ext_protate2(u,pi/2);
            y.p1 = y.p1+u*d;
            y.p2 = y.p2+u*d;     
        end
        function y = getd(obj,p)
            y = dsegment(p,[obj.p1;obj.p2]);
        end
        function y = getsd(obj,p)
            y = obj.getd(p);
            u = obj.getu;
            y = y.*sign(u(1)*p(:,2)-u(2)*p(:,1));
        end
        function y = getcta(obj,pname,side,radius)
            % getting center of tangent arc
            pname = rmspaces(pname);
            side = rmspaces(side);
            switch lower(pname)
                case 'p1'
                    u = obj.getu;
                    switch lower(side)
                        case 'l'
                            u = ext_protate2(u,pi/2);
                            y = obj.p1 + radius*u;   
                        case 'r'
                            u = ext_protate2(u,-pi/2);
                            y = obj.p1 + radius*u; 
                        otherwise
                            error('side muse be <<l>> or <<r>>.');
                    end
                case 'p2'
                    u = obj.getu;
                    switch lower(side)
                        case 'l'
                            u = ext_protate2(u,pi/2);
                            y = obj.p2 + radius*u;   
                        case 'r'
                            u = ext_protate2(u,-pi/2);
                            y = obj.p2 + radius*u; 
                        otherwise
                            error('side muse be <<l>> or <<r>>.');
                    end
                otherwise
                    error('pname must be <<p1>> or <<p2>>.');
            end  
        end
        function y = getcpa(obj,pname,side,radius)
            % get center of perpendicular arc
            pname = rmspaces(pname);
            side = rmspaces(side);
            switch lower(pname)
                case 'p1'
                    u = obj.getu;
                    switch lower(side)
                        case 'b'
                            y = obj.p1 - radius*u;   
                        case 'f'
                            y = obj.p1 + radius*u; 
                        otherwise
                            error('side muse be <<b>> or <<f>>.');
                    end
                case 'p2'
                    u = obj.getu;
                    switch lower(side)
                        case 'b'
                            y = obj.p2 - radius*u;   
                        case 'f'
                            y = obj.p2 + radius*u; 
                        otherwise
                            error('side muse be <<b>> or <<f>>.');
                    end
                otherwise
                    error('pname must be <<p1>> or <<p2>>.');
            end  
        end
        function show(obj,label,scale)
            if nargin<3
                scale = 1;
            end
            p = [obj.p1;obj.p2];
            plot(p(:,1),p(:,2));
            pc = sum(p)/2;
            vec = scale*obj.getu;
            vec = ext_protate2(vec,pi/2);
            quiver(pc(1),pc(2),vec(1),vec(2),'color','r');
            text(pc(1)+vec(1)/2,pc(2)+vec(2)/2,label);
        end
    end
end