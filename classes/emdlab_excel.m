% this mfile is developed for reading of excel files by MATLAB
% https://ComProgExpert.com

classdef emdlab_excel < handle
    
    properties
        
        data
        vars
        
    end
    
    methods
        
        function obj = emdlab_excel()
            obj.vars = struct;
        end
        
        function importExcelFile(obj, varargin)
            [~, ~, obj.data] = xlsread(varargin{:});
        end
        
        function readVariables(obj, varargin)
            
            % number of variables
            Nvars = numel(varargin);
            
            % loop over excel data
            for i = 1:size(obj.data,1)
                
                for j = 1:size(obj.data,2)
                    
                    for k = 1:Nvars
                        
                        if strcmp(varargin{k}, obj.data{i,j})
                            
                            obj.vars.(varargin{k}).value = obj.data{i,j+1};
                            
                            if isnan(obj.data{i,j+2})
                                obj.vars.(varargin{k}).unit = '';
                            else
                                obj.vars.(varargin{k}).unit = obj.data{i,j+2};
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        function str = getValueUnitString(obj, varName)
            str = [num2str(obj.vars.(varName).value), obj.vars.(varName).unit];
        end
        
        function defineAllVariablesAsGlobal(obj, em)
            
            varNames = fieldnames(obj.vars);
            for i = 1:numel(varNames)
                em.defineGlobalVariable(varNames{i}, obj.getValueUnitString(varNames{i}));
            end
              
        end
        
        function changeValueOfAllGlobalVariables(obj, em)
            
            varNames = fieldnames(obj.vars);
            for i = 1:numel(varNames)
                em.changeGlobalVariableValue(varNames{i}, obj.getValueUnitString(varNames{i}));
            end
              
        end
        
        function y = getValue(obj, varName)
            y = obj.vars.(varName).value;
        end
        
    end
    
end