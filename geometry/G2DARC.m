% arc properties class
classdef G2DARC
    properties (Constant)
        geps = 1e-6;
    end
    properties
        center
        p1
        p2
        direction
        maxLength
        maxDegree
        Nnodes
    end
    methods
        function obj = G2DARC(center,p1,p2,varargin)
            obj.center = center;
            obj.p1 = p1;
            obj.p2 = p2;
            
            obj.checkArc;
            
            p = inputParser();
            addParameter(p,'maxLength',inf);
            addParameter(p,'maxDegree',5);
            addParameter(p,'direction',1);
            addParameter(p,'Nnodes',2);
            p.CaseSensitive = 0;
            parse(p,varargin{:});
            
            obj.maxLength = p.Results.maxLength;
            obj.direction= p.Results.direction;
            obj.maxDegree = p.Results.maxDegree;
            obj.Nnodes = p.Results.Nnodes;
            
        end
        function checkArc(obj)
            % check for consistency
            p1c = obj.p1 - obj.center;
            p2c = obj.p2 - obj.center;
            r1 = norm(p1c);
            r2 = norm(p2c);
            if abs(r2-r1) > obj.geps 
                error('These points do not form an arc.');    
            end
        end
        function y = getRadius(obj)
            y = norm(obj.p1-obj.center);
        end
        function y = getAngle(obj)
            p1c = obj.p1 - obj.center;
            p2c = obj.p2 - obj.center;
            p2p1 = obj.p2 - obj.p1;
            y = 2*asin(norm(p2p1)/2/obj.getRadius);
            if obj.direction >= 0
                if det([p1c;p2c]) < 0
                    y = 2*pi - y ;
                end
            else
                if det([p1c;p2c]) > 0
                    y = 2*pi - y ;
                end
            end
        end
        function y = getip(obj)
            thetaA = obj.getAngle;
            if obj.Nnodes > 2
                n = obj.Nnodes;
            else
                % thetaP is equivalent polygonal step angle
                thetaeps = 2*asin(obj.geps/2);
                thetaP = 2*asin(obj.maxLength/2/obj.getRadius);
                n = ceil(thetaA/thetaP) + 1;
                MaxDegree = obj.maxDegree*pi/180;
                if thetaP > MaxDegree
                    tmp = floor((thetaA/MaxDegree)/thetaeps)*thetaeps;
                    n = ceil(tmp) + 1;
                end
            end
            % output points
            y = zeros(n-2,2);
            if obj.direction>=0
                stepAngle = obj.getAngle/(n-1);
            else
                stepAngle = -obj.getAngle/(n-1);
            end
            for i = 1:n-2
                y(i,:) = ext_protate2(obj.p1,i*stepAngle,obj.center);
            end
        end
        function [e,v] = getev(obj)
            v = [obj.p1;obj.getip;obj.p2];
            Nv = size(v,1);
            e = [1:Nv-1;2:Nv]';
        end
        function y = getd(obj,p)
            y = min(dsegment(p,[obj.p1;obj.getip;obj.p2]),[],2);
        end
        function y = getzd(obj,p)
            cp = [p(:,1)-obj.center(1),p(:,2)-obj.center(2)];
            index1 = abs(sqrt(sum(cp.^2,2))-obj.getRadius)<obj.geps;
%             p = p(index,:);
            u12 = obj.getu12;
            u12 = [u12(2),-u12(1)];
            index2 = p*u12'-obj.p1*u12';
            if obj.direction > 0
                index2 = index2>=-obj.geps;
            else
                index2 = index2<=obj.geps;
            end
            y = bitand(index1,index2);
        end
            
            
        function y = getmirror(obj,varargin)
            y = obj;
            y.p1 = pmirror(y.p1,varargin{:});
            y.p2 = pmirror(y.p2,varargin{:});
            y.center = pmirror(y.center,varargin{:});
            y.direction = -y.direction;
        end
        function y = getrotate(obj,varargin)
            y = obj;
            y.p1 = ext_protate2(y.p1,varargin{:});
            y.p2 = ext_protate2(y.p2,varargin{:});
            y.center = ext_protate2(y.center,varargin{:});
        end
        function y = getshift(obj,varargin)
            y = obj;
            y.p1 = pshift(y.p1,varargin{:});
            y.p2 = pshift(y.p2,varargin{:});
            y.center = pshift(y.center,varargin{:});
        end
        function y = getu1(obj)
            y = (obj.p1-obj.center)/norm(obj.p1-obj.center);
        end
        function y = getu12(obj)
            y = (obj.p2-obj.p1)/norm(obj.p2-obj.p1);
        end
        function y = getu2(obj)
            y = (obj.p2-obj.center)/norm(obj.p2-obj.center);
        end
        function y = getindent(obj,d)
            y = obj;
            u1 = y.getu1;
            u2 = y.getu2;
            y.p1 = y.p1 + d*u1;
            y.p2 = y.p2 + d*u2;
        end
        function y = getcta(obj,pname,side,radius)
            pname = rmspaces(pname);
            side = rmspaces(side);
            switch lower(pname)
                case 'p1'
                    u = obj.getu1;
                    switch lower(side)
                        case 'i'
                            y = obj.p1 - radius*u;   
                        case 'o'
                            y = obj.p1 + radius*u; 
                        otherwise
                            error('side muse be <<l>> or <<r>>.');
                    end
                case 'p2'
                    u = obj.getu2;
                    switch lower(side)
                        case 'i'
                            y = obj.p2 - radius*u;   
                        case 'o'
                            y = obj.p2 + radius*u; 
                        otherwise
                            error('side muse be <<i>> or <<o>>.');
                    end
                otherwise
                    error('pname must be <<p1>> or <<p2>>.');
            end  
        end
        function y = getcpa(obj,pname,side,radius)
            pname = rmspaces(pname);
            side = rmspaces(side);
            switch lower(pname)
                case 'p1'
                    u = obj.getu1;
                    u = ext_protate2(u,pi/2);
                    switch lower(side)
                        case 'b'    
                            y = obj.p1 - radius*u;   
                        case 'f'
                            y = obj.p1 + radius*u; 
                        otherwise
                            error('side muse be <<b>> or <<f>>.');
                    end
                case 'p2'
                    u = obj.getu2;
                    u = ext_protate2(u,pi/2);
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
            p = [obj.p1;obj.getip;obj.p2];
            plot(p(:,1),p(:,2));
            n = size(p,1);
            switch rem(n,2)
                case 0
                    xp1 = p(n/2,:);
                    xp2 = p(n/2+1,:);
                    pc = (xp1 +  xp2)/2;
                case 1
                    xp1 = p(ceil(n/2)-1,:);
                    xp2 = p(ceil(n/2)+1,:);
                    pc = p(ceil(n/2),:);
            end
            vec = scale*(xp2-xp1)/norm(xp2-xp1);
            vec = ext_protate2(vec,pi/2);
            quiver(pc(1),pc(2),vec(1),vec(2),'color','r');
            text(pc(1)+vec(1)/2,pc(2)+vec(2)/2,label);
        end
    end
end