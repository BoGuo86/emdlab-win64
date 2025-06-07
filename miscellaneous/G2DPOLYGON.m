classdef G2DPOLYGON
    properties
        v
    end
    methods
        function obj = G2DPOLYGON(v)
            obj.v = v;
        end
        function y = getd(obj,p)
            v1 = obj.v;
            v2 = circshift(v1,-1);
            y = getd(G2DLINE(v1(1,:),v2(1,:)),p);
            
            for i = 2:size(v1,1)
                y = min([y,getd(G2DLINE(v1(i,:),v2(i,:)),p)],[],2);
            end
        end
        function show(obj)
            tmp = 1:size(obj.v,1);
            tmp = tmp';
            patch('faces',[tmp,circshift(tmp,-1)],'vertices',obj.v,...
                'edgeColor','k');
        end
        function getim(obj)
            t = delaunay(obj.v);
            c = (obj.v(t(:,1),:)+obj.v(t(:,2),:)+obj.v(t(:,3),:))/3;
            d = obj.getd(c);
            t = t(d<0,:);
            
            triplot(triangulation(t,obj.v));
        end
    end
    
end