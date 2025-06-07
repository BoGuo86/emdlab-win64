% this mfile is developed for automatic generation of vbscript for generation of the design in MagNet software
% https://ComProgExpert.com

classdef emdlab_magnet < handle
    
    properties
        
        projectPath
        projectName
        designName
        fid = -1;
        
    end
    
    
    methods
        
        function obj = emdlab_magnet()
        end
        
        function delete(obj)
            
            if obj.fid > 0
                fclose(obj.fid);
            end
            
        end
        
        function openScriptFile(obj)
            
            obj.fid = fopen('MagNetScript.vbs', 'w');
            if obj.fid < 0
                MException('', 'Can not open script file.');
            end
            
            obj.initialize;
            
        end
        
        function closeAndExecuteScriptFile(obj)
            
            % close script file & run
            if fclose(obj.fid)==0
                obj.fid = -1;
            else
                MException('', 'can not close script file.');
            end
            
            magnet = actxserver('MagNet.Application');
            magnet.set('Visible',1);
            magnet.invoke('RunScript', [cd, '\MagNetScript.vbs']);

        end
        
        function initialize(obj)
            
            obj.writeLine('Call newDocument()');
            
        end
        
        function writeLine(obj, str)
            
            fprintf(obj.fid, '%s\n', str);
            
        end
        
        function newLine(obj, x1, y1, x2, y2)
            
            obj.writeLine(sprintf('Call getDocument().getView().newLine(%f, %f, %f, %f)', x1, y1, x2, y2))
            
        end
        
        function newPolyline(obj, x, y)
            
            for i = 1:length(x)-1
                obj.newLine(x(i), y(i), x(i+1), y(i+1));
            end
            
        end
        
        function newArc(obj, x1, y1, x2, y2, x3, y3)
            
            obj.writeLine(sprintf('Call getDocument().getView().newArc(%f, %f, %f, %f, %f, %f)', x1, y1, x2, y2, x3, y3))
            
        end
        
        function newCircle(obj, x, y, r)
            
            obj.writeLine(sprintf('Call getDocument().getView().newCircle(%f, %f, %f)', x, y, r))
            
        end
        
        function segmentAllEdges(obj, r)
            
            obj.writeLine(sprintf('Call getDocument().getView().selectIn(%f, %f, %f, %f, infoSetSelection, Array(infoSliceLine, infoSliceArc))', -r, -r, r, r))
            obj.writeLine('Call getDocument().getView().segmentEdges()');
            
        end
        
        function clearSketch(obj, r)
            
            obj.writeLine(sprintf('Call getDocument().getView().selectIn(%f, %f, %f, %f, infoSetSelection, Array(infoSliceLine, infoSliceArc))', -r, -r, r, r))
            obj.writeLine('Call getDocument().getView().deleteSelection()');
            
        end
        
        function selectEdge(obj, x, y)
            
            obj.writeLine(sprintf('Call getDocument().getView().selectAt(%f, %f, infoSetSelection, Array(infoSliceLine, infoSliceArc))', x, y))
            
        end
        
        function deleteSelection(obj)
            
            obj.writeLine('Call getDocument().getView().deleteSelection()');
        end
        
        function setDefaultLengthUnit(obj, str)
            
            obj.writeLine(sprintf('Call getDocument().setDefaultLengthUnit("%s")', str));
            
        end
        
        function selectAt(obj, x ,y)
            
            obj.writeLine(sprintf('Call getDocument().getView().selectAt(%f, %f, infoSetSelection)', x, y));
            
        end
        
        function makeComponentInALine(obj, x, y, depth, name, material)
            
            obj.selectAt(x,y);
            obj.writeLine(sprintf('Call getDocument().getView().makeComponentInALine(%f, Array("%s"), "NAME=%s", infoMakeComponentUnionSurfaces Or infoMakeComponentRemoveVertices)', depth, name, material));
            
        end
                
    end
    
    methods (Static)
        
        function y = cellNames2Array(names)
            y = 'Array(';
            for i = 1:numel(names)
                y = strcat(y, ['"', names{i},'"']);
            end
            y = strcat(y, ')');
        end
        
        function y = cellNames2CSV(names)
            y = '';
            for i = 1:numel(names)
                y = strcat(y, [names{i}, ',']);
            end
            y = y(1:end-1);
        end
        
        function y = getDuplicatedNamesCSV(name, i1, i2)
            
            y = name;
            for i = i1:i2
                y = strcat(y, [',', name, '_', num2str(i)]);
            end
            
        end
        
        function y = getVectorFormat(x)
            
            y = '[';
            for ii = 1:length(x)
                y = strcat(y, sprintf('%.5f,', x(ii)));
            end
            y(end) = ']';
            
        end
        
        function [xNew,yNew] = rotatePoint(x, y, a, xo, yo)
            
            if nargin == 3
                xo = 0;
                yo = 0;
            end
            
            x = x - xo;
            y = y - yo;
            
            xNew = xo + x*cos(a) - y*sin(a);
            yNew = yo + x*sin(a) + y*cos(a);
            
        end
        
        function [xNew,yNew] = rotatePointH(x, y, h, xo, yo)
            
            if nargin == 3
                xo = 0;
                yo = 0;
            end
            
            a = atan(h/norm([x-xo, y-yo]));
            
            [xNew,yNew] = emdlab_magnet.rotatePoint(x, y, a, xo, yo);
            
        end
        
    end
    
end