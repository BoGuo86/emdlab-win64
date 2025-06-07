% this mfile is developed for automatic generation of vbscript for generation of the design in ANSYS Maxwell software
% https://ComProgExpert.com

classdef emdlab_maxwell < handle
    
    properties
        
        projectPath
        projectName
        designName
        fid
        isFileOpen (1,1) logical = false;
        
    end
    
    methods
        
        function obj = emdlab_maxwell()
        end
        
        function delete(obj)
            
            if obj.isFileOpen
                fclose(obj.fid);
            end
            obj.isFileOpen = false;
            
        end
        
        function openScriptFile(obj, index)
            
            % index: 0 for Maxwell, 1 for ANSYS Electronics, default value is zero
            
            if nargin == 1
                index = 0;
            end
            
            obj.fid = fopen('MaxwellScript.vbs', 'w');
            if obj.fid < 0
                error('Can not open file.');
            end
            
            obj.isFileOpen = true;
            obj.initialize(index);
            
        end
        
        function closeAndExecuteScriptFile(obj)
            
            % close script file & run
            fclose(obj.fid);
            obj.isFileOpen = false;
            system('MaxwellScript.vbs');
            
        end
        
        function writeLine(obj, str)
            
            fprintf(obj.fid, '%s\n', str);
            
        end
        
        function defineGlobalVariable(obj, varName, varValue)
            
            fprintf(obj.fid, 'oProject.ChangeProperty Array("NAME:AllTabs", Array("NAME:ProjectVariableTab", Array("NAME:PropServers",  _\n');
            fprintf(obj.fid, '"ProjectVariables"), Array("NAME:NewProps", Array("NAME:$%s", "PropType:=", "VariableProp", "UserDef:=",  _\n', varName);
            fprintf(obj.fid, 'true, "Value:=", "%s"))))\n', varValue);
            
        end
        
        function changeGlobalVariableValue(obj, varName, newValue)
            
            fprintf(obj.fid, 'oProject.ChangeProperty Array("NAME:AllTabs", Array("NAME:ProjectVariableTab", Array("NAME:PropServers",  _\n');
            fprintf(obj.fid, '"ProjectVariables"), Array("NAME:ChangedProps", Array("NAME:$%s", "Value:=", "%s"))))\n', varName, newValue);
            
        end
        
        function defineLocalVariable(obj, varName, varValue)
            
            fprintf(obj.fid, 'oDesign.ChangeProperty Array("NAME:AllTabs", Array("NAME:LocalVariableTab", Array("NAME:PropServers",  _\n');
            fprintf(obj.fid, '"LocalVariables"), Array("NAME:NewProps", Array("NAME:%s", "PropType:=", "VariableProp", "UserDef:=",  _\n', varName);
            fprintf(obj.fid, 'true, "Value:=", "%s"))))\n', varValue);
            
        end
        
        function newProject(obj)
            
            obj.writeLine('Set oProject = oDesktop.NewProject');
            
        end
        
        function insertDesign(obj, type, name, solutionType)
            
            obj.writeLine(sprintf('oProject.InsertDesign "%s", "%s", "%s", ""', type, name, solutionType));
            obj.writeLine(sprintf('Set oDesign = oProject.SetActiveDesign("%s")', name));
            
        end
        
        function drawCircle(obj, name, xc, yc, r)
            
            xc = emdlab_maxwell.change2str(xc);
            yc = emdlab_maxwell.change2str(yc);
            r = emdlab_maxwell.change2str(r);
            
            obj.set3DModelerAsActiveEditor;
            obj.writeLine('oEditor.CreateCircle Array("NAME:CircleParameters", "IsCovered:=", true, "XCenter:=",  _');
            obj.writeLine(sprintf('"%s", "YCenter:=", "%s", "ZCenter:=", "0mm", "Radius:=", "%s", "WhichAxis:=",  _', xc, yc, r));
            obj.writeLine(sprintf('"Z", "NumSegments:=", "0"), Array("NAME:Attributes", "Name:=", "%s", "Flags:=",  _', name));
            obj.writeLine('"", "Color:=", "(132 132 193)", "Transparency:=", 0, "PartCoordinateSystem:=",  _');
            obj.writeLine('"Global", "UDMId:=", "", "MaterialValue:=", "" & Chr(34) & "vacuum" & Chr(34) & "", "SolveInside:=",  _');
            obj.writeLine('  true)');
            
        end
        
        function drawRectangle(obj, name, x0, y0, w, h)
            
            x0 = emdlab_maxwell.change2str(x0);
            y0 = emdlab_maxwell.change2str(y0);
            w = emdlab_maxwell.change2str(w);
            h = emdlab_maxwell.change2str(h);
            
            obj.set3DModelerAsActiveEditor;
            obj.writeLine('oEditor.CreateRectangle Array("NAME:RectangleParameters", "IsCovered:=", true, "XStart:=",  _');
            obj.writeLine(sprintf('"%s", "YStart:=", "%s", "ZStart:=", "0mm", "Width:=", "%s", "Height:=",  _', x0, y0, w));
            obj.writeLine(sprintf('"%s", "WhichAxis:=", "Z"), Array("NAME:Attributes", "Name:=", "%s", "Flags:=",  _', h, name));
            obj.writeLine('"", "Color:=", "(132 132 193)", "Transparency:=", 0, "PartCoordinateSystem:=",  _');
            obj.writeLine('"Global", "UDMId:=", "", "MaterialValue:=", "" & Chr(34) & "vacuum" & Chr(34) & "", "SolveInside:=",  true)');
            
        end
        
        function subtract(obj, blankParts, toolParts, keepOriginals)
            
            if nargin == 3, keepOriginals = 'false'; end
            keepOriginals = emdlab_maxwell.change2str(keepOriginals);
            
            obj.set3DModelerAsActiveEditor;
            obj.writeLine(sprintf('oEditor.Subtract Array("NAME:Selections", "Blank Parts:=", "%s", "Tool Parts:=",  _', blankParts));
            obj.writeLine(sprintf('"%s"), Array("NAME:SubtractParameters", "KeepOriginals:=", %s)', toolParts, keepOriginals));
            
        end
        
        function intersect(obj, names, keepOriginals)
            
            if nargin == 2, keepOriginals = 'false'; end
            keepOriginals = emdlab_maxwell.change2str(keepOriginals);
            
            obj.set3DModelerAsActiveEditor;
            obj.writeLine(sprintf('oEditor.Intersect Array("NAME:Selections", "Selections:=", "%s"), Array("NAME:IntersectParameters", "KeepOriginals:=", %s)', ...
                names, keepOriginals));
            
        end
        
        function openPolyline(obj)
            
            obj.set3DModelerAsActiveEditor;
            fprintf(obj.fid, 'oEditor.CreatePolyline Array("NAME:PolylineParameters", "IsPolylineCovered:=", true, "IsPolylineClosed:=", true, _\n');
            fprintf(obj.fid, 'Array("NAME:PolylinePoints" _\n');
            
        end
        
        function addPoint2Polyline(obj, x, y)
            
            x = emdlab_maxwell.change2str(x);
            y = emdlab_maxwell.change2str(y);
            
            fprintf(obj.fid, ', Array("NAME:PLPoint", "X:=", "%s", "Y:=", "%s", "Z:=", "0mm") _\n', x, y);
            
        end
        
        function startPolylineSegments(obj)
            
            fprintf(obj.fid, '), Array("NAME:PolylineSegments" _\n');
            
        end
        
        function addLineSegment2Polyline(obj, startIndex)
            
            fprintf(obj.fid, ', Array("NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", %s, "NoOfPoints:=", 2) _\n', startIndex);
            
        end
        
        function addArcSegment2Polyline(obj, startIndex, n)
            
            fprintf(obj.fid, ', Array("NAME:PLSegment", "SegmentType:=", "Arc", "StartIndex:=", %s, "NoOfPoints:=", 3, "NoOfSegments:=", "%s") _\n', startIndex, n);
            
        end
        
        function closePolyline(obj, name)
            
            fprintf(obj.fid, '), Array("NAME:PolylineXSection", "XSectionType:=", "None", "XSectionOrient:=",  _\n');
            fprintf(obj.fid, '"Auto", "XSectionWidth:=", "0mm", "XSectionTopWidth:=", "0mm", "XSectionHeight:=",  _\n');
            fprintf(obj.fid, '"0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=", "Corner")), Array("NAME:Attributes", "Name:=",  _\n');
            fprintf(obj.fid, '"%s", "Flags:=", "", "Color:=", "(132 132 193)", "Transparency:=", 0, "PartCoordinateSystem:=",  _\n', name);
            fprintf(obj.fid, '"Global", "UDMId:=", "", "MaterialValue:=", "" & Chr(34) & "vacuum" & Chr(34) & "", "SolveInside:=", true)\n');
            
        end
        
        function duplicateMirror(obj, names, xb, yb, zb, xn, yn, zn)
            
            xb = emdlab_maxwell.change2str(xb);
            yb = emdlab_maxwell.change2str(yb);
            zb = emdlab_maxwell.change2str(zb);
            xn = emdlab_maxwell.change2str(xn);
            yn = emdlab_maxwell.change2str(yn);
            zn = emdlab_maxwell.change2str(zn);
            
            obj.set3DModelerAsActiveEditor;
            fprintf(obj.fid, 'oEditor.DuplicateMirror Array("NAME:Selections", "Selections:=", "%s", "NewPartsModelFlag:=",  _\n', names);
            fprintf(obj.fid, '"Model"), Array("NAME:DuplicateToMirrorParameters", "DuplicateMirrorBaseX:=", "%s", "DuplicateMirrorBaseY:=",  _\n', xb);
            fprintf(obj.fid, '"%s", "DuplicateMirrorBaseZ:=", "%s", "DuplicateMirrorNormalX:=", "%s", "DuplicateMirrorNormalY:=",  _\n', yb, zb, xn);
            fprintf(obj.fid, '"%s", "DuplicateMirrorNormalZ:=", "%s"), Array("NAME:Options", "DuplicateAssignments:=", false)\n', yn, zn);
            
        end
        
        function mirror(obj, names, xb, yb, zb, xn, yn, zn)
            
            xb = emdlab_maxwell.change2str(xb);
            yb = emdlab_maxwell.change2str(yb);
            zb = emdlab_maxwell.change2str(zb);
            xn = emdlab_maxwell.change2str(xn);
            yn = emdlab_maxwell.change2str(yn);
            zn = emdlab_maxwell.change2str(zn);
            
            obj.set3DModelerAsActiveEditor;
            fprintf(obj.fid, 'oEditor.Mirror Array("NAME:Selections", "Selections:=", "%s", "NewPartsModelFlag:=",  _\n', names);
            fprintf(obj.fid, '"Model"), Array("NAME:MirrorParameters", "MirrorBaseX:=", "%s", "MirrorBaseY:=",  _\n', xb);
            fprintf(obj.fid, '"%s", "MirrorBaseZ:=", "%s", "MirrorNormalX:=", "%s", "MirrorNormalY:=",  _\n', yb, zb, xn);
            fprintf(obj.fid, '"%s", "MirrorNormalZ:=", "%s")\n', yn, zn);
            
        end
        
        function rotate(obj, name, angle)
            
            obj.set3DModelerAsActiveEditor;
            fprintf(obj.fid, 'oEditor.Rotate Array("NAME:Selections", "Selections:=", "%s", "NewPartsModelFlag:=",  _\n', name);
            fprintf(obj.fid, '"Model"), Array("NAME:RotateParameters", "RotateAxis:=", "Z", "RotateAngle:=",  _\n');
            fprintf(obj.fid, '"%s")\n', angle);
            
        end
        
        function unite(obj, names, keepOriginals)
            
            if nargin == 2, keepOriginals = 'false'; end
            keepOriginals = emdlab_maxwell.change2str(keepOriginals);
            
            obj.set3DModelerAsActiveEditor;
            fprintf(obj.fid, 'oEditor.Unite Array("NAME:Selections", "Selections:=", "%s"), Array("NAME:UniteParameters", "KeepOriginals:=", %s)\n', ...
                names, keepOriginals);
            
        end
        
        function duplicateAroundAxis(obj, names, axisName, d, n)
            
            obj.set3DModelerAsActiveEditor;
            fprintf(obj.fid, 'oEditor.DuplicateAroundAxis Array("NAME:Selections", "Selections:=", "%s", "NewPartsModelFlag:=",  _\n', names);
            fprintf(obj.fid, '"Model"), Array("NAME:DuplicateAroundAxisParameters", "CreateNewObjects:=", true, "WhichAxis:=",  _\n');
            fprintf(obj.fid, '"%s", "AngleStr:=", "%s", "NumClones:=", "%s"), Array("NAME:Options", "DuplicateAssignments:=", false)\n', axisName, d, n );
            
        end
        
        function changeObjectName(obj, oldName, newName)
            
            obj.set3DModelerAsActiveEditor;
            fprintf(obj.fid, 'oEditor.ChangeProperty Array("NAME:AllTabs", Array("NAME:Geometry3DAttributeTab", Array("NAME:PropServers",  _\n');
            fprintf(obj.fid, '"%s"), Array("NAME:ChangedProps", Array("NAME:Name", "Value:=", "%s"))))\n', oldName, newName);
            
        end
        
        function changeObjectColor(obj, name, r, g, b)
            
            obj.set3DModelerAsActiveEditor;
            fprintf(obj.fid, 'oEditor.ChangeProperty Array("NAME:AllTabs", Array("NAME:Geometry3DAttributeTab", Array("NAME:PropServers",  _\n');
            fprintf(obj.fid, '"%s"), Array("NAME:ChangedProps", Array("NAME:Color", "R:=", %d, "G:=", %d, "B:=", %d))))\n', name, r, g, b);
            
        end
        
        function changeObjectTransparent(obj, name, value)
            
            obj.set3DModelerAsActiveEditor;
            fprintf(obj.fid, 'oEditor.ChangeProperty Array("NAME:AllTabs", Array("NAME:Geometry3DAttributeTab", Array("NAME:PropServers",  _\n');
            fprintf(obj.fid, '"%s"), Array("NAME:ChangedProps", Array("NAME:Transparent", "Value:=", %f))))\n', name, value);
            
        end
        
        function changeObjectMaterial(obj, name, materialName)
            
            obj.set3DModelerAsActiveEditor;
            obj.writeLine('oEditor.ChangeProperty Array("NAME:AllTabs", Array("NAME:Geometry3DAttributeTab", Array("NAME:PropServers",  _');
            obj.writeLine(sprintf('"%s"), Array("NAME:ChangedProps", Array("NAME:Material", "Value:=", "" & Chr(34) & "%s" & Chr(34) & ""))))', name, materialName));
            
        end
        
        function assignCurrent(obj, name, value, direction)
            
            obj.writeLine('Set oModule = oDesign.GetModule("BoundarySetup")');
            obj.writeLine(sprintf('oModule.AssignCurrent Array("NAME:%s", "Objects:=", Array("%s"), "Current:=",  _', name, name));
            obj.writeLine(sprintf('"%s", "IsPositive:=", %s)\n', value, direction));
            
        end
        
        function createEquationCurve2D01(obj, name, xt, yt)
            
            obj.writeLine('oEditor.CreateEquationCurve Array("NAME:EquationBasedCurveParameters", "XtFunction:=",  _');
            obj.writeLine(sprintf('"%s", "YtFunction:=", "%s", "ZtFunction:=", "0", "tStart:=", "0", "tEnd:=",  _', xt, yt));
            obj.writeLine('"1", "NumOfPointsOnCurve:=", "0", "Version:=", 1, Array("NAME:PolylineXSection", "XSectionType:=",  _');
            obj.writeLine('"None", "XSectionOrient:=", "Auto", "XSectionWidth:=", "0", "XSectionTopWidth:=",  _');
            obj.writeLine('"0", "XSectionHeight:=", "0", "XSectionNumSegments:=", "0", "XSectionBendType:=",  _');
            obj.writeLine(sprintf('"Corner")), Array("NAME:Attributes", "Name:=", "%s", "Flags:=", "", "Color:=",  _', name));
            obj.writeLine('"(132 132 193)", "Transparency:=", 0, "PartCoordinateSystem:=", "Global", "UDMId:=",  _');
            obj.writeLine('"", "MaterialValue:=", "" & Chr(34) & "vacuum" & Chr(34) & "", "SolveInside:=",  _');
            obj.writeLine('true)');
            
        end
        
        function drawSegmentByEC(obj, name, x0, y0, x1, y1)
            
            xt = [x0, '+ _t * (', x1, '-', x0, ')'];
            yt = [y0, '+ _t * (', y1, '-', y0, ')'];
            
            obj.createEquationCurve2D01(name, xt, yt);
            
        end
        
        function assignTorque(obj, names)
            
            obj.writeLine('Set oModule = oDesign.GetModule("MaxwellParameterSetup")');
            names = emdlab_maxwell.csv2array(names);
            obj.writeLine(sprintf('oModule.AssignTorque Array("NAME:Static_Torque", "Coordinate System:=", "Global", "Is Positive:=", true, "Objects:=", %s)', names));    
            
        end
        
    end
    
    methods (Access = private)
        
        function initialize(obj, index)
            
            obj.writeLine('Dim oAnsoftApp');
            obj.writeLine('Dim oDesktop');
            obj.writeLine('Dim oProject');
            obj.writeLine('Dim oDesign');
            obj.writeLine('Dim oEditor');
            obj.writeLine('Dim oModule');
            
            if index == 0
                obj.writeLine('Set oAnsoftApp = CreateObject("AnsoftMaxwell.MaxwellScriptInterface")');
            else
                obj.writeLine('Set oAnsoftApp = CreateObject("Ansoft.ElectronicsDesktop")');
            end
            
            obj.writeLine('Set oDesktop = oAnsoftApp.GetAppDesktop()');
            obj.writeLine('oDesktop.RestoreWindow');
            
        end
        
        function set3DModelerAsActiveEditor(obj)
            
            obj.writeLine('Set oEditor = oDesign.SetActiveEditor("3D Modeler")');
            
        end
        
    end
    
    methods (Static)
        
        function y = cellNames2Array(names)
            
            y = 'Array(';
            for i = 1:numel(names)-1
                y = strcat(y, ['"', names{i},'",']);
            end
            y = strcat(y, ['"', names{end},'")']);
            
        end
        
        function y = csv2array(x)
            
            x = strsplit(x, ',');
            y = emdlab_maxwell.cellNames2Array(x);
            
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
        
        function y = change2str(x)
            
            if isnumeric(x)
                
                y = num2str(x);
                
            elseif islogical(x)
                
                if x
                    y = 'true';
                else
                    y = 'false';
                end
                
            else
                
                y = x;
                
            end
            
        end
        
    end
    
end