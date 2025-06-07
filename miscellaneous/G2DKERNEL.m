classdef G2DKERNEL < handle 
    properties (Constant = true)
        % geometry epsilon for length
        gleps = 1e-6 
        % geometry epsilon for angle
        gaeps = 2*atan(1e-6/2) 
    end  
    properties (SetAccess = private)
        % key points
        Vertices
        % edges: segments, arcs and ...
        Edges
        % closed boundaries
        Loops
        % 2D regions
        Faces 
    end
    methods    
        function obj = G2DKERNEL() 
            obj.Vertices = struct();
            obj.Edges = struct();
            obj.Loops = struct();
            obj.Faces = struct();      
        end  
        %% Key Point Functions
        function newkp(obj,kpname,kpc)
            % new key point
            % kpname: key point name
            % kpc: key point coordinate
            if ~ischar(kpname)
                error('kpname must be a string.');
            end
            kpname = obj.rmspaces(kpname);
            [m,n] = size(kpc); 
            if ~(m==1&&n==2)
                error('kpc must be a 2x1 matrix.');
            end
            if isfield(obj.Vertices,kpname)
                error('another kp with the same name exist');
            end
            obj.Vertices.(kpname) = kpc;
        end
        function dkpname = newdkp(obj,kpc)
            % new default key point
            % kpc: key point coordinate
            kpnames = fieldnames(obj.Vertices);
            for i = 1:numel(kpnames)
                if norm(obj.Vertices.(kpnames{i})-kpc)<obj.gleps
                    dkpname = kpnames{i};
                    return
                end
            end
            dkpname = obj.getdkpname();
            obj.newkp(dkpname,kpc);
        end    
        function dkpname = getdkpname(obj)
            % get default key point name
            index = 1;
            while isfield(obj.Vertices,(['kp',num2str(index)]))
                index = index + 1;
            end        
            dkpname = ['kp',num2str(index)];            
        end  
        function kpnames = newdkps(obj,p)
            % new default key points
            kpnames = cell(1,size(p,1));
            for i = 1:size(p,1)
                [obj,kpnames{i}] = obj.newdkp(p(i,:));
            end     
        end        
        %% Edge Functions 
        % Segment Functions
        function news(obj,sname,kp1name,kp2name,varargin)
            % removing extra spaces
            sname = obj.rmspaces(sname);
            kp1name = obj.rmspaces(kp1name);
            kp2name = obj.rmspaces(kp2name);
            % defenition of new line
            obj.Edges.(sname).type = 'segment';
            obj.Edges.(sname).kp1name = kp1name;
            obj.Edges.(sname).kp2name = kp2name;
            obj.Edges.(sname).value = G2DSEGMENT(obj.Vertices.(kp1name),...
                obj.Vertices.(kp2name),varargin{:});
        end
        function dsname = newds(obj,kp1name,kp2name,varargin)
            % removing extra spaces
            dsname = obj.getdsname;
            kp1name = obj.rmspaces(kp1name);
            kp2name = obj.rmspaces(kp2name);
            % defenition of new line
            obj.Edges.(dsname).type = 'segment';
            obj.Edges.(dsname).kp1name = kp1name;
            obj.Edges.(dsname).kp2name = kp2name;
            obj.Edges.(dsname).value = G2DSEGMENT(obj.Vertices.(kp1name),...
                obj.Vertices.(kp2name),varargin{:});
        end
        function dsname = getdsname(obj)
            % get default segment name
            index = 1;
            while isfield(obj.Edges,(['S',num2str(index)]))
                index = index + 1;
            end
            dsname = ['S',num2str(index)];        
        end 
        function [dsname,dkp1name,dkp2name] = newdswdkps(obj,p1,p2,varargin)
            % new default line with default key points
            dkp1name = obj.newdkp(p1);
            dkp2name = obj.newdkp(p2);
            dsname = obj.newds(dkp1name,dkp2name,varargin{:});
        end
        % Arc Functions
        function newarccpp(obj,aname,kp1name,kp2name,center,varargin)
            % getting point coordinates
            aname = obj.rmspaces(aname);
            kp1name = obj.rmspaces(kp1name);
            kp2name = obj.rmspaces(kp2name);
            % defenition of new arc
            obj.Edges.(aname).type = 'arc';
            obj.Edges.(aname).kp1name = kp1name;
            obj.Edges.(aname).kp2name = kp2name;
            obj.Edges.(aname).value = G2DARC(center,...
                    obj.Vertices.(kp1name),obj.Vertices.(kp2name),varargin{:});
        end
        function darcname = getdaname(obj)
            % get default arc name
            sindex = 1;
            while isfield(obj.Edges,(['A',num2str(sindex)]))
                sindex = sindex + 1;
            end 
            darcname = ['A',num2str(sindex)];        
        end
        function darcname = newdarccpp(obj,varargin)
            % getting deefault arc name
            darcname = obj.getdaname();
            % defenition of arc
            obj.newarccpp(darcname,varargin{:});
        end      
        function [daname,dkp1name,dkp2name] = newdarccppwdkps(obj,center,p1,p2,...
                varargin)
            % defenition of default kps
            dkp1name = obj.newdkp(p1);
            dkp2name = obj.newdkp(p2);
            % getting default arc name
            daname = obj.getdaname();
            % defenition of arc
            obj.newarccpp(daname,dkp1name,dkp2name,center,varargin{:});
        end  
        % Transform Functions
        function cmirrore(obj,nename,ename,varargin)
            % removing extra spaces
            nename = obj.rmspaces(nename);
            ename = obj.rmspaces(ename);
            % checking for existence of tool edge
            if isfield(obj.Edges,nename)
                error(['Another edge with the name <<',nename,'>> exist'])
            else
                obj.Edges.(nename) = obj.Edges.(ename);
                obj.Edges.(nename).value = ...
                    obj.Edges.(nename).value.getmirror(varargin{:});
                obj.Edges.(nename).kp1name = obj.newdkp(obj.Edges.(nename).value.p1);
                obj.Edges.(nename).kp2name = obj.newdkp(obj.Edges.(nename).value.p2);
            end
        end   
        function cshifte(obj,nename,ename,varargin)
            % removing extra spaces
            nename = obj.rmspaces(nename);
            ename = obj.rmspaces(ename);
            % checking for existence of tool edge
            if isfield(obj.Edges,nename)
                error(['Another edge with the name <<',nename,'>> exist'])
            else
                obj.Edges.(nename) = obj.Edges.(ename);
                obj.Edges.(nename).value = ...
                    obj.Edges.(nename).value.getshift(varargin{:});
                obj.Edges.(nename).kp1name = obj.newdkp(obj.Edges.(nename).value.p1);
                obj.Edges.(nename).kp2name = obj.newdkp(obj.Edges.(nename).value.p2);
            end
        end 
        function crotatee(obj,nename,ename,varargin)
            % removing extra spaces
            nename = obj.rmspaces(nename);
            ename = obj.rmspaces(ename);
            % checking for existence of tool boundary
            if isfield(obj.Edges,nename)
                error(['Another boundary with the name <<',nename,'>> exist'])
            else
                obj.Edges.(nename) = obj.Edges.(ename);
                obj.Edges.(nename).value = ...
                    obj.Edges.(nename).value.getrotate(varargin{:});
                obj.Edges.(nename).kp1name = obj.newdkp(obj.Edges.(nename).value.p1);
                obj.Edges.(nename).kp2name = obj.newdkp(obj.Edges.(nename).value.p2);
            end
        end        
        % Copy and Transform Functions
        function nbname = cmirrored(obj,ename,varargin)
            % removing extra spaces
            ename = obj.rmspaces(ename);
            switch obj.Edges.(ename).type
                case 'segment'
                    nbname = obj.getdsname;
                case 'arc'
                    nbname = obj.getdaname;
            end
            obj.cmirrore(nbname,ename,varargin{:});
        end
        function nbname = cshifted(obj,ename,varargin)
            % removing extra spaces
            ename = obj.rmspaces(ename);
            switch obj.Edges.(ename).type
                case 'segment'
                    nbname = obj.getdsname;
                case 'arc'
                    nbname = obj.getdaname;
            end
            obj.cshifte(nbname,ename,varargin{:});
        end
        function nbname = crotateed(obj,ename,varargin)
            ename = obj.rmspaces(ename);
            switch obj.Edges.(ename).type
                case 'segment'
                    nbname = obj.getdsname;
                case 'arc'
                    nbname = obj.getdaname;
            end
            obj.crotatee(nbname,ename,varargin{:});
        end
        % Arc Segment Joining Functions
        function [obj,darcname] = newdarctl(obj,sname,pname,side,radius,angle,direction,varargin)
            % removing extra spaces
            sname = obj.rmspaces(sname);
            pname = obj.rmspaces(pname);
            % checking for existance of line
            if ~isfield(obj.Edges,sname)
                error(['Segment with the name <<',sname,'>> does not exist.']);
            end
            % getting center of tangent arc
            center = obj.Edges.(sname).value.getcta(pname,side,radius);
            switch lower(pname)
                case 'p1'
                    p1arc = obj.Edges.(sname).value.p1;
                    p2arc = ext_protate2(p1arc,sign(direction)*angle*pi/180,center);
                    [obj,darcname] = obj.newdarccppwdkps(center,p1arc,p2arc,...
                        'direction',direction,varargin{:});
                case 'p2'
                    p1arc = obj.Edges.(sname).value.p2;
                    p2arc = ext_protate2(p1arc,sign(direction)*angle*pi/180,center);
                    [obj,darcname] = obj.newdarccppwdkps(center,p1arc,p2arc,...
                        'direction',direction,varargin{:});
                otherwise
                    error('pname must be <<p1>> or <<p2>>.');
            end
        end
        function [obj,darcname] = newdarcta(obj,lname,pname,side,radius,angle,direction,varargin)
            lname = obj.rmspaces(lname);
            if ~ismember(lname,fieldnames(obj.Edges))
                error(['line with the name <<',lname,'>> does not exist.']);
            end
            center = obj.Edges.(lname).value.getcta(pname,side,radius);
            switch lower(pname)
                case 'p1'
                    p1arc = obj.Edges.(lname).value.p1;
                    p2arc = ext_protate2(p1arc,sign(direction)*angle*pi/180,center);
                    [obj,darcname] = obj.newdarccppwdkps(center,p1arc,p2arc,direction,varargin{:});
                case 'p2'
                    p1arc = obj.Edges.(lname).value.p2;
                    p2arc = ext_protate2(p1arc,sign(direction)*angle*pi/180,center);
                    [obj,darcname] = obj.newdarccppwdkps(center,p1arc,p2arc,direction,varargin{:});
                otherwise
                    error('pname must be <<p1>> or <<p2>>.');
            end
        end
        function [obj,darcname] = newdarcpl(obj,lname,pname,side,radius,angle,direction,varargin)
            lname = obj.rmspaces(lname);
            if ~ismember(lname,fieldnames(obj.Edges))
                error(['line with the name <<',lname,'>> does not exist.']);
            end
            center = obj.Edges.(lname).value.getcpa(pname,side,radius);
            switch lower(pname)
                case 'p1'
                    p1arc = obj.Edges.(lname).value.p1;
                    p2arc = ext_protate2(p1arc,sign(direction)*angle*pi/180,center);
                    [obj,darcname] = obj.newdarccppwdkps(center,p1arc,p2arc,'direction',direction,varargin{:});
                case 'p2'
                    p1arc = obj.Edges.(lname).value.p2;
                    p2arc = ext_protate2(p1arc,sign(direction)*angle*pi/180,center);
                    [obj,darcname] = obj.newdarccppwdkps(center,p1arc,p2arc,'direction',direction,varargin{:});
                otherwise
                    error('pname must be <<p1>> or <<p2>>.');
            end
        end
        function [obj,darcname] = newdarcpa(obj,lname,pname,side,radius,angle,direction,varargin)
            lname = obj.rmspaces(lname);
            if ~ismember(lname,fieldnames(obj.Edges))
                error(['line with the name <<',lname,'>> does not exist.']);
            end
            center = obj.Edges.(lname).value.getcpa(pname,side,radius);
            switch lower(pname)
                case 'p1'
                    p1arc = obj.Edges.(lname).value.p1;
                    p2arc = ext_protate2(p1arc,sign(direction)*angle*pi/180,center);
                    [obj,darcname] = obj.newdarccppwdkps(center,p1arc,p2arc,direction,varargin{:});
                case 'p2'
                    p1arc = obj.Edges.(lname).value.p2;
                    p2arc = ext_protate2(p1arc,sign(direction)*angle*pi/180,center);
                    [obj,darcname] = obj.newdarccppwdkps(center,p1arc,p2arc,direction,varargin{:});
                otherwise
                    error('pname must be <<p1>> or <<p2>>.');
            end
        end
        function [obj,nl1name,nl2name] = splitline(obj,lname)
            lname = obj.rmspaces(lname);
            if ~ismember(lname,fieldnames(obj.Edges))
                error(['line with the name <<',lname,'>> does not exist.']);
            end
            p1 = obj.Edges.(lname).value.p1;
            p2 = obj.Edges.(lname).value.p2;
            maxLength = obj.Edges.(lname).value.maxLength;
            obj.Edges = rmfield(obj.Edges,lname);
            [obj,nl1name] = obj.newdswdkps(p1,(p1+p2)/2,maxLength);
            [obj,nl2name] = obj.newdswdkps((p1+p2)/2,p2,maxLength);
        end
        %% Loop Functions
        function dlbname = getdloopname(obj)
            % get default loop name
            lindex = 1;
            while isfield(obj.Loops,(['L',num2str(lindex)]))
                lindex = lindex + 1;
            end       
            dlbname = ['L',num2str(lindex)];          
        end  
        function newloop(obj,cbname,varargin)
            % new closed boundary
            % bs: boundaries
            % ds: directions
            cbname = obj.rmspaces(cbname);
            nbs = numel(varargin);
            if rem(nbs,2) == 1
                error('Number of boundaries and directions must be even.');
            end
            obj.Loops.(cbname).bs = varargin(1:2:end);
            obj.Loops.(cbname).ds = varargin(2:2:end);
        end    
        function y = getloopbps(obj,lname)
            lname = obj.rmspaces(lname);
            y = zeros([],2);
            for i = 1:numel(obj.Loops.(lname).bs)
                b = obj.Edges.(obj.rmspaces(obj.Loops.(lname).bs{i})).value;
                d = obj.Loops.(lname).ds{i};
                p = [b.p1;b.getip;b.p2];
                if d<0
                    p = flipud(p);
                end
                y = [y;p(1:end-1,:)];
            end
        end 
        function [E,V] = getloopev(obj,lname)
            lname = obj.rmspaces(lname);
            Nledges = numel(obj.Loops.(lname).bs);
            e = cell(1,Nledges);
            Ne = zeros(1,Nledges);
            v = cell(1,Nledges);
            Nv = zeros(1,Nledges);
            for i = 1:Nledges
                [e{i},v{i}] = obj.Edges.(obj.rmspaces(obj.Loops.(lname).bs{i})).value.getev;
                Nv(i) = size(v{i},1);
                Ne(i) = Nv(i)-1;
                if obj.Loops.(lname).ds{i} < 0
                    e{i} = fliplr(e{i});
                end
            end
            E = zeros(sum(Ne),2);
            V = zeros(sum(Nv),2);
            Eindex = 0;
            Vindex = 0;
            for i = 1:Nledges
                E(Eindex+1:Eindex+Ne(i),:) = e{i}+Vindex;
                Eindex = Eindex + Ne(i);
                V(Vindex+1:Vindex+Nv(i),:) = v{i};
                Vindex = Vindex + Nv(i);
            end
        end
        function y = getloopbpsm(obj,lname)
            lname = obj.rmspaces(lname);
            y = zeros([],2);
            for i = 1:numel(obj.Loops.(lname).bs)
                b = obj.Edges.(obj.rmspaces(obj.Loops.(lname).bs{i})).value;
                d = obj.Loops.(lname).ds{i};
                if strcmpi(obj.Edges.(obj.rmspaces(obj.Loops.(lname).bs{i})).type,'segment')
                    p = [b.p1;b.p2];
                else
                    p = [b.p1;b.getip;b.p2];
                end
                if d<0
                    p = preverse(p);
                end
                y = [y;p(1:end-1,:)];
            end
        end 
        function lbname = newdloop(obj,varargin)
            % new defualt name closed bound
            lbname = obj.getdloopname;           
            obj.newloop(lbname,varargin{:});            
        end      
        function newcloop(obj,lname,center,radius,varargin)
            % new circular closed bound
            % varargin = (maxLength,maxDegree)
            % defenition of dkps
            dkp1name = obj.newdkp(center+[radius,0]);
            dkp2name = obj.newdkp(center-[radius,0]);
            % defenition of darcs
            a1name = obj.newdarccpp(dkp1name,dkp2name,center,varargin{:});
            a2name = obj.newdarccpp(dkp2name,dkp1name,center,varargin{:});
            % defenition of closed bound
            obj.newloop(lname,a1name,1,a2name,1);
        end   
        function [obj,dcbname] = newdcircularcb(obj,c,r,varargin)   
            dcbname = obj.getdloopname;
            obj.newcircularcb(dcbname,c,r,varargin{:});
        end         
        function obj = newploop(obj,lname,p)
            % new polygonal loop
            [obj,kpnames] = obj.newdkps(p);
            % boundaries and directions
            bsds = cell(1,2*numel(kpnames));          
            for i = 1:(numel(kpnames)-1)
                [obj,bsds{2*i-1}] = obj.newds(kpnames{i},kpnames{i+1});
                bsds{2*i} = 1;
            end
            [obj,bsds{end-1}] = obj.newds(kpnames{end},kpnames{1});
            bsds{end} = 1;     
            obj = obj.newloop(lname,bsds{:});           
        end        
        function [obj,dcbname] = newdpolygonalcb(obj,p,MaxLength)
            dcbname = obj.getdloopname;
            obj = obj.newploop(dcbname,p,MaxLength);
        end           
        function obj = newloopwd(obj,cbname,h0,varargin)
            % new closed bound with domain
            obj = obj.newloop(cbname,varargin{:});
            obj = obj.newd(cbname,h0,cbname);        
        end       
        %% Face Functions
        function newface(obj,fname,varargin)
            % new face
            if isempty(varargin)
                error('At lest one boundary must be specefied.')
            end
            % removing extra spaces
            fname = obj.rmspaces(fname);
            % defenition of new domain
            obj.Faces.(fname) = varargin;
        end
        function [f,v] = getfacefvv(obj,fname)
            % removing extra spaces
            fname = obj.rmspaces(fname);
            % getting Xbps and Ybps of domain
            Nloops = numel(obj.Faces.(fname));
            Xbps = cell(1,Nloops);
            Ybps = cell(1,Nloops);
            Nv = 0;
            for i = 1:Nloops
                p = obj.getloopbps(obj.Faces.(fname){i});
                Nv = Nv + size(p,1);
                if i == 1
                    [x,y] = poly2cw(p(:,1),p(:,2));
                else
                    [x,y] = poly2ccw(p(:,1),p(:,2));
                end
                Xbps{i} = x;
                Ybps{i} = y;
            end
            % calculation of f and v of domain
            % generation of inititial initial mesh
            % constraint edges
            ces = zeros(Nv,2);
            temp = 0;
            for i = 1:numel(Xbps)
                Np = size(Xbps{i},1);
                ces(temp+1:temp+Np,:) = [temp+1:temp+Np;[temp+2:temp+Np,temp+1]]';
                temp = temp + Np;
            end
            [xbps,ybps] = polyjoin(Xbps,Ybps);
            pfix = [xbps(~isnan(xbps)),ybps(~isnan(xbps))];
            % making triangulation
            dt = delaunayTriangulation(pfix,ces);
            % removing outer triangles
            c = dt.incenter;
            sd = dpoly(c,[Xbps{1},Ybps{1}]);
            for i = 2:numel(Xbps)
                sd = max([sd,-dpoly(c,[Xbps{i},Ybps{i}])],[],2);
            end
            % adding face to database
            f = dt.ConnectivityList(sd<-obj.gleps,:);
            v = dt.Points;
        end
        function [f,v] = getfacefv(obj,fname)
            % removing extra spaces
            fname = obj.rmspaces(fname);
            [E,V] = getfaceev(obj,fname);
            % making triangulation
            dt = delaunayTriangulation(V,E);
            % removing outer triangles
            c = dt.incenter;
            sd = ext_dpoly2d(c,E,V);
            % adding face to database
            f = dt.ConnectivityList(sd<-obj.gleps,:);
            v = dt.Points;
        end
        function [E,V] = getfaceev(obj,fname)
            % removing extra spaces
            fname = obj.rmspaces(fname);
            % getting Xbps and Ybps of domain
            Nloops = numel(obj.Faces.(fname));
            e = cell(1,Nloops);
            Ne = zeros(1,Nloops);
            v = cell(1,Nloops);
            Nv = zeros(1,Nloops);
            for i = 1:Nloops
                [e{i},v{i}] = obj.getloopev(obj.Faces.(fname){i});
                Ne(i) = size(e{i},1);
                Nv(i) = size(v{i},1);
            end
            E = zeros(sum(Ne),2);
            V = zeros(sum(Nv),2);
            Eindex = 0;
            Vindex = 0;
            for i = 1:1
                E(Eindex+1:Eindex+Ne(i),:) = e{i}+Vindex;
                Eindex = Eindex + Ne(i);
                V(Vindex+1:Vindex+Nv(i),:) = v{i};
                Vindex = Vindex + Nv(i);
            end
            for i = 2:Nloops
                E(Eindex+1:Eindex+Ne(i),:) = fliplr(e{i})+Vindex;
                Eindex = Eindex + Ne(i);
                V(Vindex+1:Vindex+Nv(i),:) = v{i};
                Vindex = Vindex + Nv(i);
            end
            % unification
            [V,~,ic] = uniquetol(V,obj.gleps,'ByRows',true);
            E = ic(E);
        end
        function [f,v] = getfacefvm(obj,fname)
            % removing extra spaces
            fname = obj.rmspaces(fname);
            % getting Xbps and Ybps of domain
            Nloops = numel(obj.Faces.(fname));
            Xbps = cell(1,Nloops);
            Ybps = cell(1,Nloops);
            Nv = 0;
            for i = 1:Nloops
                p = obj.getloopbpsm(obj.Faces.(fname){i});
                Nv = Nv + size(p,1);
                if i == 1
                    [x,y] = poly2cw(p(:,1),p(:,2));
                else
                    [x,y] = poly2ccw(p(:,1),p(:,2));
                end
                Xbps{i} = x;
                Ybps{i} = y;
            end
            % calculation of f and v of domain
            % generation of inititial initial mesh
            % constraint edges
            ces = zeros(Nv,2);
            temp = 0;
            for i = 1:numel(Xbps)
                Np = size(Xbps{i},1);
                ces(temp+1:temp+Np,:) = [temp+1:temp+Np;[temp+2:temp+Np,temp+1]]';
                temp = temp + Np;
            end
            [xbps,ybps] = polyjoin(Xbps,Ybps);
            pfix = [xbps(~isnan(xbps)),ybps(~isnan(xbps))];
            % making triangulation
            dt = delaunayTriangulation(pfix,ces);
            % removing outer triangles
            c = dt.incenter;
            sd = dpoly(c,[Xbps{1},Ybps{1}]);
            for i = 2:numel(Xbps)
                sd = max([sd,-dpoly(c,[Xbps{i},Ybps{i}])],[],2);
            end
            % adding face to database
            f = dt.ConnectivityList(sd<-obj.gleps,:);
            v = dt.Points;
        end
        function mz = getdmmz(obj,fname)
            % removing extra spaces
            fname = obj.rmspaces(fname);
            % edge vertices constraint
            [e,v] = obj.getfaceev(fname);         
            % calling mesh generator
            mz = MeshGenerator0(v,e);
            % mesh smoothing
            mz = mz.moveNodes;     
        end
        function newcface(obj,fname,c,r,varargin)
            % new circular domain
            obj.newcloop(fname,c,r,varargin{:});
            obj.newface(fname,fname);  
        end     
        function newpface(obj,fname,p)  
            % new polygonal face
            obj.newploop(fname,p);
            obj.newface(fname,fname);
        end      
        function newrectangulard(obj,dname,h0,x,y,w,h,hw,hh)
            % define needed key points
            kpnames = obj.newdkps([x,y;x+w,y;x+w,y+h;x,y+h]);
            L1 = obj.newds(kpnames{1},kpnames{2},hw);
            L2 = obj.newds(kpnames{2},kpnames{3},hh);
            L3 = obj.newds(kpnames{3},kpnames{4},hw);
            L4 = obj.newds(kpnames{4},kpnames{1},hh);
            obj.newloop(dname,{L1,L2,L3,L4},[1 1 1 1]);
            obj.newd(dname,h0,dname);  
        end              
        %% Tools For Modifying Geometry  
        % Fillets
        function daname = filletss(obj,s1name,s2name,r,varargin)
            % removing extra spaces
            s1name = obj.rmspaces(s1name);
            s2name = obj.rmspaces(s2name);
            % checking for existance of lines
            if ~(isfield(obj.Edges,s1name)&& ... 
                    strcmpi(obj.Edges.(s1name).type,'segment'))
                error(['Segment with name <<',s1name,'>> does not exist'])
            end
            if ~(isfield(obj.Edges,s2name)&& ... 
                    strcmpi(obj.Edges.(s2name).type,'segment'))
                error(['Segment with name <<',s2name,'>> does not exist'])
            end
            % getting name of line kps
            s1kps = {obj.Edges.(s1name).kp1name,obj.Edges.(s1name).kp2name};
            s2kps = {obj.Edges.(s2name).kp1name,obj.Edges.(s2name).kp2name};
            % finding comman kp
            [jkps,i1,i2] = intersect(s1kps,s2kps);
            % checking for jointing line at one kp
            if numel(jkps) == 2
                error('Segments are coincident')
            elseif numel(jkps) == 0
                error('Segments do not intersect at comman key point')
            else
                l1 = obj.Edges.(s1name).value.getindent(r);
                l2 = obj.Edges.(s2name).value.getindent(r);
                c = obj.getssi(l1,l2);
                if c == false
                    l1 = obj.Edges.(s1name).value.getindent(-r);
                    l2 = obj.Edges.(s2name).value.getindent(r);
                    c = obj.getssi(l1,l2);
                end
                if c == false
                    l1 = obj.Edges.(s1name).value.getindent(r);
                    l2 = obj.Edges.(s2name).value.getindent(-r);
                    c = obj.getssi(l1,l2);
                end
                if c == false
                    l1 = obj.Edges.(s1name).value.getindent(-r);
                    l2 = obj.Edges.(s2name).value.getindent(-r);
                    c = obj.getssi(l1,l2);
                end
                if c == false
                    error('Fillet is not possible.');
                end 
            end
            p1 = obj.getPPOS(obj.Edges.(s1name).value,c);
            p2 = obj.getPPOS(obj.Edges.(s2name).value,c);
            % new key points
            kp1name = obj.newdkp(p1);
            kp2name = obj.newdkp(p2); 
            % modifying lines
            % line 1
            if i1 == 1
                obj.Edges.(s1name).kp1name = kp1name;
                obj.Edges.(s1name).value.p1 = p1;
            else
                obj.Edges.(s1name).kp2name = kp1name;
                obj.Edges.(s1name).value.p2 = p1;
            end
            % line 2
            if i2 == 1
                obj.Edges.(s2name).kp1name = kp2name;
                obj.Edges.(s2name).value.p1 = p2;
            else
                obj.Edges.(s2name).kp2name = kp2name;
                obj.Edges.(s2name).value.p2 = p2;
            end
            % adding new arc
            u1 = p1 - obj.Vertices.(jkps{1});
            u2 = p2 - obj.Vertices.(jkps{1});
            if det([u1;u2])<0
                daname = obj.newdarccppwdkps(c,p1,p2,varargin{:});
            elseif det([u1;u2])>0
                daname = obj.newdarccppwdkps(c,p2,p1,varargin{:});
            end
        end    
        function [obj,daname] = filletsa(obj,sname,aname,r,varargin)
            % removing extra spaces
            sname = obj.rmspaces(sname);
            aname = obj.rmspaces(aname);
            % checking for existance of lines
            if ~(isfield(obj.Edges,sname)&&strcmpi(obj.Edges.(sname).type,'segment'))
                error(['line with name <<',sname,'>> do not exist.'])
            end
            if ~(isfield(obj.Edges,aname)&&strcmpi(obj.Edges.(aname).type,'arc'))
                error(['arc with name <<',aname,'>> do not exist.'])
            end
            % getting name of line kps
            l1kps = {obj.Edges.(sname).kp1name,obj.Edges.(sname).kp2name};
            l2kps = {obj.Edges.(aname).kp1name,obj.Edges.(aname).kp2name};
            % finding comman kp
            [jkps,i1,i2] = intersect(l1kps,l2kps);
            % checking for jointing line at one kp
            if numel(jkps) == 2
                error('segment and arc are coincident at to point')
            elseif numel(jkps) == 0
                error('lines do not intersect at comman key point')
            else
                s = obj.Edges.(sname).value.getindent(r);
                a = obj.Edges.(aname).value.getindent(r);
                p = obj.getsai(s,a);
                if p == false
                    s = obj.Edges.(sname).value.getindent(-r);
                    a = obj.Edges.(aname).value.getindent(r);
                    p = obj.getsai(s,a);
                end
                if p == false
                    s = obj.Edges.(sname).value.getindent(r);
                    a = obj.Edges.(aname).value.getindent(-r);
                    p = obj.getsai(s,a);
                end
                if p == false
                    s = obj.Edges.(sname).value.getindent(-r);
                    a = obj.Edges.(aname).value.getindent(-r);
                    p = obj.getsai(s,a);
                end
                if p == false
                    error('Fillet is not possible.');
                end
                    
                p1 = obj.getPPOS(obj.Edges.(sname).value,p);
                p2 = p-a.center;
                p2 = obj.getsai(GLINEPC(a.center,p+a.getRadius*p2),obj.Edges.(aname).value);
                
                
                
                % new key points
                [obj,kp1name] = obj.newdkp(p1);
                [obj,kp2name] = obj.newdkp(p2);
                % modifying line and arc
                if i1 == 1
                    obj.Edges.(sname).kp1name = kp1name;
                    obj.Edges.(sname).value.p1 = p1;
                else
                    obj.Edges.(sname).kp2name = kp1name;
                    obj.Edges.(sname).value.p2 = p1;
                end
                
                 if i2 == 1
                    obj.Edges.(aname).kp1name = kp2name;
                    obj.Edges.(aname).value.p1 = p2;
                else
                    obj.Edges.(aname).kp2name = kp2name;
                    obj.Edges.(aname).value.p2 = p2;
                 end
                
                 u1 = p1-obj.Vertices.(jkps{1});
                 u2 = p2-obj.Vertices.(jkps{1});
                
                if det([u1;u2])<0
                    [obj,daname] = obj.newdarccppwdkps(p,p1,p2,varargin{:});
                  elseif det([u1;u2])>0
                    [obj,daname] = obj.newdarccppwdkps(p,p2,p1,varargin{:});
                     end      
            end   
        end 
        function y = findRayCircleIntersect(obj,p,u,c,r)
            d = norm(p-c);
            u = u/norm(u);
            cp = p-c;
            if d<r-obj.gleps
                c1 = 1;
                c2 = 2*cp*u';
                c3 = norm(cp)^2-r^2;
                beta1 = (-c2+sqrt(c2^2-4*c1*c3))/2/c1;  
                beta2 = (-c2-sqrt(c2^2-4*c1*c3))/2/c1; 
                y = p + max(beta1,beta2)*u;
            elseif abs(d-r)<obj.gleps
                if (cp*u')>=0
                    y = p;
                end
            else
                error('Not impelemented yet.');
            end
        end
        function y = getAAI(obj,a1,a2)
            u = a2.center - a1.center;
            v = [u(2),-u(1)];
            r1 = a1.getRadius;
            r2 = a2.getRadius;
            
            r12 = norm(u);
            
            delta = -(r12^2-(r1+r2)^2)*(r12^2-(r1-r2)^2);
            t = sqrt(delta/4/r12^2);
            s = 0.5*((r1^2-r2^2)/r12^2+1);
            t = sqrt(r1^2/r12^2-s^2);
            
            if r12<obj.gleps
                y = false;
                return
            end
            
            if delta<0
                y = false;
                return
            end
            
            if abs(r12-r1-r2)<obj.gleps
                y = a1.center+(r1/(r1+r2))*u;
                if ~(a1.getzd(y)&&a2.getzd(y))
                    y = false;
                end
            elseif abs(r12-abs(r1-r2))<obj.gleps
                y = a1.center+(r1/(r1-r2))*u;
                if ~(a1.getzd(y)&&a2.getzd(y))
                    y = false;
                end
            else
                y = [a1.center + s*u + t*v;a1.center + s*u - t*v];
                i1 = a1.getzd(y);
                i2 = a2.getzd(y);
                y = y(bitand(i1,i2),:);
                if isempty(y)
                    y = false;
                end
            end
            
        end
        function y = getSAI(obj,l,a)
            
            Delta = l.p1-a.center;
            D = l.p2 - l.p1;
            delta = (D*Delta')^2 - norm(D)^2*(norm(Delta).^2-a.getRadius^2);
            
            if delta<-obj.gleps
                y = false;
                return
            end
            
            t1 = (-(D*Delta')+sqrt(delta))/norm(D)^2;
            t2 = (-(D*Delta')-sqrt(delta))/norm(D)^2;
            
            if t1 >=0 && t1 <=1
                y1 = l.p1+t1*(l.p2-l.p1);
                s = det([y1-a.p1;a.getu12]);
                 if a.direction >0
                     if s<0
                         y1 = false;
                     end
                 else
                    if s>0
                         y1 = false;
                    end
                end
            else
                y1 = false;
            end
            
            if t2 >=0 && t2 <=1
                y2 = l.p1+t2*(l.p2-l.p1);
                s = det([y2-a.p1;a.getu12]);
                 if a.direction >0
                     if s<0
                         y2 = false;
                     end
                 else
                    if s>0
                         y2 = false;
                    end
                end
            else
                y2 = false;
            end
            
            if y1 == false
                if y2 == false
                    y = false;
                else
                    y = y2;
                end
            elseif y2 == false
                y = y1;
            else
                y = false;
            end
            
        end
        function y = getPPOS(obj,s,p)
            % point projection on segment
            if ~isa(s,'G2DSEGMENT')
                error('Input must be a segment.')
            end
            a = s.p1;
            b = s.p2;
            c = p;
            u = s.getu;
            u = ext_protate2(u,pi/2);
            d = p+u;
            
            A = [-b(1)+a(1) d(1)-c(1)
                -b(2)+a(2) d(2)-c(2)];
            bb = [a(1)-c(1);a(2)-c(2)];

            y = A\bb;
            if y(1)>1 || y(1) <0
                y = false;
                return;
            end
            
            y = a+y(1)*(b-a);
            
        end   
        function [obj,daname] = filletaa(obj,a1name,a2name,r,varargin)
            % removing extra spaces
            a1name = obj.rmspaces(a1name);
            a2name = obj.rmspaces(a2name);
            % checking for existance of lines
            if ~(isfield(obj.Edges,a1name)&&strcmpi(obj.Edges.(a1name).type,'arc'))
                error(['arc with name <<',a1name,'>> do not exist.'])
            end
            if ~(isfield(obj.Edges,a2name)&&strcmpi(obj.Edges.(a2name).type,'arc'))
                error(['arc with name <<',a2name,'>> do not exist.'])
            end
            % getting name of line kps
            l1kps = {obj.Edges.(a1name).kp1name,obj.Edges.(a1name).kp2name};
            l2kps = {obj.Edges.(a2name).kp1name,obj.Edges.(a2name).kp2name};
            % finding comman kp
            [jkps,i1,i2] = intersect(l1kps,l2kps);
            % checking for jointing line at one kp
            if numel(jkps) == 2
                error('lines are coincident')
            elseif numel(jkps) == 0
                error('lines do not intersect at comman key point')
            else
                a1 = obj.Edges.(a1name).value.getindent(r);
                a2 = obj.Edges.(a2name).value.getindent(r);
                p = obj.getAAI(a1,a2);
                if p == false
                    a1 = obj.Edges.(a1name).value.getindent(-r);
                    a2 = obj.Edges.(a2name).value.getindent(r);
                    p = obj.getAAI(a1,a2);
                end
                if p == false
                    a1 = obj.Edges.(a1name).value.getindent(r);
                    a2 = obj.Edges.(a2name).value.getindent(-r);
                    p = obj.getAAI(a1,a2);
                end
                if p == false
                    a1 = obj.Edges.(a1name).value.getindent(-r);
                    a2 = obj.Edges.(a2name).value.getindent(-r);
                    p = obj.getAAI(a1,a2);
                end
                if p == false
                    error('Fillet is not possible.');
                end
                    
                p1 = p-a1.center;
                p1 = obj.getsai(GLINEPC(a1.center,p+a1.getRadius*p1),...
                    obj.Edges.(a1name).value);
                
                if p1 == false
                    error('Fillet is not possible.');
                end
                
                p2 = p-a2.center;
                p2 = obj.getsai(GLINEPC(a2.center,p+a2.getRadius*p2),...
                    obj.Edges.(a2name).value);
                
                if p2 == false
                    error('Fillet is not possible.');
                end
                
                
                
                % new key points
                [obj,kp1name] = obj.newdkp(p1);
                [obj,kp2name] = obj.newdkp(p2);
                % modifying line and arc
                if i1 == 1
                    obj.Edges.(a1name).kp1name = kp1name;
                    obj.Edges.(a1name).value.p1 = p1;
                else
                    obj.Edges.(a1name).kp2name = kp1name;
                    obj.Edges.(a1name).value.p2 = p1;
                end
                
                 if i2 == 1
                    obj.Edges.(a2name).kp1name = kp2name;
                    obj.Edges.(a2name).value.p1 = p2;
                else
                    obj.Edges.(a2name).kp2name = kp2name;
                    obj.Edges.(a2name).value.p2 = p2;
                 end
                
                 u1 = p1-obj.Vertices.(jkps{1});
                 u2 = p2-obj.Vertices.(jkps{1});
                
                if det([u1;u2])<0
                    [obj,daname] = obj.newdarccppwdkps(p,p1,p2,varargin{:});
                    %                 obj = obj.news(l1name,p1name{1},kp1name,'maxLength',h1);
                    %                 obj = obj.news(l2name,kp2name,p2name{1},'maxLength',h2);
                    %                 [obj,daname] = obj.newdarccpp(kp1name,kp2name,pc+u*y,1,varargin{:});
                elseif det([u1;u2])>0
                    [obj,daname] = obj.newdarccppwdkps(p,p2,p1,varargin{:});
                    
                    %                 obj = obj.news(l1name,kp1name,p1name{1},'maxLength',h1);
                    %                 obj = obj.news(l2name,p2name{1},kp2name,'maxLength',h2);
                    %                 [obj,daname] = obj.newdarccpp(kp1name,kp2name,pc+u*y,1,varargin{:});
                end
                
                
            end
            
        end
        % Trim Functions
        function [obj,sname] = trimss(obj,s1name,s2name,Dir)
            if nargin<4
                Dir = 1;
            end
            % removing extra spaces
            s1name = obj.rmspaces(s1name);
            s2name = obj.rmspaces(s2name);
            % checking for existance of lines
            if ~(isfield(obj.Edges,s1name)&& ... 
                    strcmpi(obj.Edges.(s1name).type,'segment'))
                error(['Segment with name <<',s1name,'>> does not exist'])
            end
            if ~(isfield(obj.Edges,s2name)&& ... 
                    strcmpi(obj.Edges.(s2name).type,'segment'))
                error(['Segment with name <<',s2name,'>> does not exist'])
            end
            s1 = obj.Edges.(s1name).value;
            s2 = obj.Edges.(s2name).value;
            % get ssi
            c = obj.getssi(s1,s2);
            vec1 = s1.getu;
            vec2 = s2.p1 - s1.p1;
            % trim
            if c ~= false
                obj.Edges = rmfield(obj.Edges,s2name);

                if det([vec1;vec2])>0
                    if Dir>0
                        [obj,sname] = obj.newdswdkps(c,s2.p2);
                    else
                        [obj,sname] = obj.newdswdkps(s2.p1,c);
                    end
                else
                    if Dir>0
                        [obj,sname] = obj.newdswdkps(s2.p1,c);
                    else
                        [obj,sname] = obj.newdswdkps(c,s2.p2);
                    end
                end
            end
        end
        % Area and Length
        function y = getFaceArea(obj,fname)
            fname = obj.rmspaces(fname);
            if obj.checkFaceName(fname)
                [f,v] = obj.getfacefvm(fname);
                vec2 = v(f(:,2),:)-v(f(:,1),:);
                v13 = v(f(:,3),:)-v(f(:,1),:);
                y = 0.5*(vec2(:,1).*v13(:,2)-vec2(:,2).*v13(:,1));
                y = sum(y);
            else
                error(['Face with name <<',fname,'>> does not exist']);
            end
        end
        %% Geometry Visualization
        function plotkps(obj)
            % plot key points
            axis off equal
            hold all
            kpnames = fieldnames(obj.Vertices);
            for i = 1:numel(kpnames)
                p = obj.Vertices.(kpnames{i});
                plot(p(1),p(2),'*','color','r');
                text(p(1),p(2),kpnames{i});
            end
        end
        function plotes(obj,scale)
            % setting scale
            if nargin < 2
                scale = 1;
            end
            axis off equal
            hold all
            % getting boundary names
            bnames = fieldnames(obj.Edges);
            % plotting all boundaries
            for i = 1:numel(bnames)
                obj.Edges.(bnames{i}).value.show(bnames{i},scale);
            end
        end 
        function plotwf(obj)
            % getting boundary names
            bnames = fieldnames(obj.Edges);
            % plotting all boundaries with innerpoints
            hold all
            axis off equal
            for i = 1:numel(bnames)
                p = [obj.Edges.(bnames{i}).value.p1;...
                    obj.Edges.(bnames{i}).value.getip;...
                    obj.Edges.(bnames{i}).value.p2];
                plot(p(:,1),p(:,2));
                plot(p(:,1),p(:,2),'.','color','k');
            end
            obj.plotkps;
            zoom on
        end
        function plotfs(obj)
            % plot domains
            axis off equal
            hold all
            % getting domain names
            fnames = fieldnames(obj.Faces);  
            % plotting each domain
            for i = 1:numel(fnames)
                [f,v] = obj.getfacefv(fnames{i});
                patch('Faces', f,'Vertices',v ,...
                'FaceColor',rand(1,3),'EdgeColor', 'none','FaceAlpha',0.5);
            end  
            zoom on;
        end  
        function plotim(obj)
            % plot initial mesh
            axis off equal
            hold all
            % getting domain names
            fnames = fieldnames(obj.Faces);  
            % plotting each domain
            for i = 1:numel(fnames)
                [f,v] = obj.getfacefv(fnames{i});
                patch('Faces', f,'Vertices',v ,...
                'FaceColor',rand(1,3),'EdgeColor', 'k');
            end  
            zoom on;
        end  
        function plotimm(obj)
            % plot initial mesh
            axis off equal
            hold all
            % getting domain names
            fnames = fieldnames(obj.Faces);  
            % plotting each domain
            for i = 1:numel(fnames)
                [f,v] = obj.getfacefvm(fnames{i});
                patch('Faces', f,'Vertices',v ,...
                'FaceColor',rand(1,3),'EdgeColor', 'k');
            end  
            zoom on;
        end 
        %% Import and Export  
        function exportstl(obj,fdir)
            % getting domain names
            dnames = fieldnames(obj.Faces);
            % exporting each domain geometry
            for i = 1:numel(dnames)
                % opening file
                f = fopen([fdir,'\',dnames{i},'.STL'],'w');
                % writing domain name
                fwrite(f,dnames{i},'*char');
                for j = length(dnames{i})+1:80
                    fwrite(f,' ','*char');
                end
                % writing number of domain triangles
                fwrite(f,uint32(obj.Faces.(dnames{i}).Nf),'uint32');
                % writing each triangle points
                for j = 1:obj.Faces.(dnames{i}).Nf
                    faces = obj.Faces.(dnames{i}).f(j,:);
                    vetices = obj.Faces.(dnames{i}).v(faces,:);
                    fwrite(f,faces,'single');
                    fwrite(f,[vetices(1,:),0],'single');
                    fwrite(f,[vetices(2,:),0],'single');
                    fwrite(f,[vetices(3,:),0],'single');
                    fwrite(f,0.5,'uint16');
                end
                % closing file
                fclose(f);
            end
        end   
        %% Editing Functions
        function obj = editEdgeName(obj,newNname,oldName)
            newNname = obj.rmspaces(newNname);
            oldName = obj.rmspaces(oldName);
            obj.Edges.(newNname) = obj.Edges.(oldName);
            obj.Edges = rmfield(obj.Edges,oldName);
        end
        function y = checkFaceName(obj,fname)
            fname = obj.rmspaces(fname);
            y = isfield(obj.Faces,fname);
        end
    end 
    methods (Static = true)
        function y = getssi(s1,s2)
            % get segment segment intersection
            % check for class type of inputs
            if ~(isa(s1,'G2DSEGMENT')&&isa(s2,'G2DSEGMENT'))
                error('Inputs must be segments.')
            end
            u1 = s1.getu;
            d1 = det([u1;s2.p1-s1.p1]);
            d2 = det([u1;s2.p2-s1.p1]);
            if d1*d2>0
                y = false;
                return;
            else
                u2 = s2.getu;
                d1 = det([u2;s1.p1-s2.p1]);
                d2 = det([u2;s1.p2-s2.p1]);
                if d1*d2>0
                    y = false;
                    return;
                else
                    u2 = [u2(2);-u2(1)];
                    s = (s2.p1-s1.p1)*u2/(u1*u2);
                    y = s1.p1 + s*u1;
                end
            end
        end
        function y = getsai(s,a)
            % segment arc intersection
            Delta = s.p1-a.center;
            D = s.p2 - s.p1;
            delta = (D*Delta')^2 - norm(D)^2*(norm(Delta).^2-a.getRadius^2);
            % not intersect
            if delta<0
                y = false;
            elseif delta==0
                t = -(D*Delta')/norm(D)^2;
                if t>=0 && t<=1
                    y = s.p1+t*(s.p2-s.p1);
                else
                    y = false;
                end
            else
                t1 = (-(D*Delta')+sqrt(delta))/norm(D)^2;
                t2 = (-(D*Delta')-sqrt(delta))/norm(D)^2;
                if t1 >=0 && t1 <=1
                    y1 = s.p1+t1*(s.p2-s.p1);
                    s = det([y1-a.p1;a.getu12]);
                    if a.direction >0
                        if s<0
                            y1 = false;
                        end
                    else
                        if s>0
                            y1 = false;
                        end
                    end
                else
                    y1 = false;
                end
                if t2 >=0 && t2 <=1
                    y2 = s.p1+t2*(s.p2-s.p1);
                    s = det([y2-a.p1;a.getu12]);
                    if a.direction >0
                        if s<0
                            y2 = false;
                        end
                    else
                        if s>0
                            y2 = false;
                        end
                    end
                else
                    y2 = false;
                end
                if y1 == false
                    if y2 == false
                        y = false;
                    else
                        y = y2;
                    end
                elseif y2 == false
                    y = y1;
                else
                    y = false;
                end
            end
        end
        function y = rmspaces(x)
            % removing extra spaces in an experission
            if ~ischar(x)
                error('Input must be string.');
            end
            y = strjoin(strsplit(x,' '),'');
        end
    end
end