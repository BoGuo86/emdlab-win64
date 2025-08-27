
% data needeeed for reference elements

function y = getedata(ename)

switch upper(ename)
    
    case 'TL3'
        y.B = @(s,t) [1,s,t];
        y.dBs = [0,1,0];
        y.dBt = [0,0,1];
        y.M = [1,0,0;1,1,0;1,0,1];
        y.G = @(s,t) [1,s,t]/y.M;
        y.GG = [y.dBs;y.dBt]/y.M;
        
    case 'TL6'
        y.B = @(s,t) [1,s,t,s.^2,s.*t,t.^2];
        y.dBs = @(s,t) [0,1,0,2*s,t,0];
        y.dBt = @(s,t) [0,0,1,0,s,2*t];
        s = [0;1;0;0.5;0.5;0];
        t = [0;0;1;0;0.5;0.5];
        y.M = [ones(6,1),s,t,s.^2,s.*t,t.^2];
        y.G = @(s,t) [1,s,t,s.^2,s.*t,t.^2]/y.M;
        y.GG = @(s,t) [y.dBs(s,t);y.dBt(s,t)]/y.M;
        
    case 'QL4'
        y.B = @(s,t) [1,s,t,s.*t];
        y.dBs = @(s,t) [0,1,0,t];
        y.dBt = @(s,t) [0,0,1,s];
        s = [-1;1;1;-1];
        t = [-1;-1;1;1];
        y.M = [ones(4,1),s,t,s.*t];
        y.G = @(s,t) [1,s,t,s.*t]/y.M;
        y.GG = @(s,t) [y.dBs(s,t);y.dBt(s,t)]/y.M;
        
    case 'TTL4'
        y.B = @(u,v,w) [1,u,v,w];
        y.dBu = [0,1,0,0];
        y.dBv = [0,0,1,0];
        y.dBw = [0,0,0,1];
        u = [0;1;0;0];
        v = [0;0;1;0];
        w = [0;0;0;1];
        y.M = [ones(4,1),u,v,w];
        y.G = @(u,v,w) [1,u,v,w]/y.M;
        y.GG = [y.dBu;y.dBv;y.dBw]/y.M;
        
    otherwise
        error('this element is not in element data library');
end

end
        
