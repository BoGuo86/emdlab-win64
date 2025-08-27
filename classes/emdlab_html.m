classdef emdlab_html < handle
    
    properties
        
        fid
        fname
        
    end
    
    
    methods
        
        function obj = emdlab_html()
        end
        
        function openHTMLFile(obj, fileName)
            
            if nargin == 1
                obj.fname = 'emdlab.html';
            else
                obj.fname = [fileName, '.html'];
            end
            
            obj.fid = fopen(obj.fname, 'w');
            if obj.fid < 0
                error('Can not open file.');
            end
            
        end
        
        function closeAndExecuteHTMLFile(obj)
            
            % close script file & run
            fclose(obj.fid);
            web(obj.fname);
            
        end
        
        function writeLine(obj, varargin)
            
            fprintf(obj.fid, varargin{:});
            
        end
        
        function writeH1(obj, varargin)
            
            obj.writeLine('<h1>');
            obj.writeLine(varargin{:});
            obj.writeLine('</h1>');
            
        end
        
        function writeH2(obj, varargin)
            
            obj.writeLine('<h2>');
            obj.writeLine(varargin{:});
            obj.writeLine('</h2>');
            
        end
        
        function writeH3(obj, varargin)
            
            obj.writeLine('<h3>');
            obj.writeLine(varargin{:});
            obj.writeLine('</h3>');
            
        end
        
        function writeH4(obj, varargin)
            
            obj.writeLine('<h4>');
            obj.writeLine(varargin{:});
            obj.writeLine('</h4>');
            
        end
        
        function writeH5(obj, varargin)
            
            obj.writeLine('<h5>');
            obj.writeLine(varargin{:});
            obj.writeLine('</h5>');
            
        end
        
        function beginHTML(obj)
            
            obj.writeLine('<html>');
            obj.setStyle;
            
        end
        
        function setStyle(obj)
            
            obj.writeLine('<style>');
            obj.writeLine('table, th, td {border: 1px solid black;}');
            obj.writeLine('th, td {padding: 5px;}');
            obj.writeLine('</style>');
            
        end
        
        function endHTML(obj)
            obj.writeLine('</html>');
        end
        
        function beginBody(obj)
            obj.writeLine('<body>');
        end
        
        function endBody(obj)
            obj.writeLine('</body>');
        end
        
        function beginTable(obj)
            obj.writeLine('<table>');
        end
        
        function endTable(obj)
            obj.writeLine('</table>');
        end
        
        function writeTableHeader(obj, varargin)
            
            obj.writeLine('<tr>');
            
            for i = 1:numel(varargin)
                obj.writeLine('<th>');
                obj.writeLine(varargin{i})
                obj.writeLine('</th>');
            end
            
            obj.writeLine('</tr>');
            
        end
        
        function writeTableRow(obj, varargin)
            
            obj.writeLine('<tr>');
            
            for i = 1:numel(varargin)
                obj.writeLine('<td>');
                obj.writeLine(varargin{i})
                obj.writeLine('</td>');
            end
            
            obj.writeLine('</tr>');
            
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
        
    end
    
end