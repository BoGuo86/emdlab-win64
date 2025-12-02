% emdlab data base class for 3D geometries

classdef emdlab_g3d_db < handle

    properties

        % stp enities
        entities = {};

        % unique id of entities
        ids (1,:) double;

        cadPath = '';

    end

    properties (Dependent = true)
        % total number of entities
        Ne (1,1) double;
    end

    methods

        %% constructor and destructor
        function obj = emdlab_g3d_db()
%             % Get full path of the running file
%             filePath   = mfilename('fullpath');
%             folderPath = fileparts(filePath);
% 
%             % Define cad folder path
%             cadFolder = fullfile(folderPath, 'cad');
% 
%             % If cad folder exists, remove it
%             if exist(cadFolder, 'dir')
%                 rmdir(cadFolder, 's');  % 's' removes folder and all contents
%             end
% 
%             % Create cad folder again
%             mkdir(cadFolder);
        end

        function id = getUniqueID(obj)
        end

        function y = get.Ne(obj)
            y = numel(obj.entities);
        end

        function stpPath = getEntitySTPPathByName(obj, name)
            for i = 1:obj.Ne
                if strcmpi(name, obj.entities{i}.name)
                    stpPath = obj.entities{i}.stpPath;
                    return;
                end
            end
            error('Wrong entity name.');
        end

        function [entityID, entityHandle] = addBox(obj, name, x0, y0, z0, X, Y, Z)

            fid = fopen("C:\emdlab-win64\tmp\emdlab_gmsh_geoFile.geo", 'w'); 
            fprintf(fid,'SetFactory("OpenCASCADE");\n');            
            fprintf(fid,'Box(1) = {%f, %f, %f, %f, %f, %f};\n', x0, y0, z0, X, Y, Z);
            fclose(fid);
            obj.runGeoSaveStepSTL(name);            

        end

        function [entityID, entityHandle] = addCylinder(obj, name, x0, y0, z0, X, Y, Z, r)

            fid = fopen("C:\emdlab-win64\tmp\emdlab_gmsh_geoFile.geo", 'w'); 
            fprintf(fid,'SetFactory("OpenCASCADE");\n');            
            fprintf(fid,'Cylinder(1) = {%f, %f, %f, %f, %f, %f, %f, 2*Pi};\n', x0, y0, z0, X, Y, Z, r);
            fclose(fid);
            obj.runGeoSaveStepSTL(name);            

        end

        function extrude(obj, g, faceName, z1, z2)

            g.write_geo_file;
            index = g.getFaceIndexByTag(faceName);

            % define a new geo file           
            fid = fopen("C:\emdlab-win64\tmp\emdlab_gmsh_geoFile.geo", 'a');
            fprintf(fid, "Extrude {0, 0, %.16f} {Surface{%s};}\n", z2-z1, join(string(1:numel(g.faces)),','));            
            fprintf(fid, "Translate {0, 0, %.16f} {Volume{%d};}\n", z1, index);
            fprintf(fid, "Recursive Delete { Volume{%s};}\n", join(string(setdiff(1:numel(g.faces), index)),','));
            fprintf(fid, 'Coherence;\n');            
            fclose(fid);

            obj.runGeoSaveStepSTL(faceName); 
            
        end

        function addFace(obj, g, faceName, newFaceName)

            if nargin < 4, newFaceName = faceName; end

            g.write_geo_file;
            index = g.getFaceIndexByTag(faceName);

            % define a new geo file           
            fid = fopen("C:\emdlab-win64\tmp\emdlab_gmsh_geoFile.geo", 'a');
            fprintf(fid, "Recursive Delete { Surface{%s};}\n", join(string(setdiff(1:numel(g.faces), index)),','));
            fprintf(fid, 'Coherence;\n');            
            fclose(fid);

            obj.runGeoSaveStepSTL(newFaceName); 
            
        end

        function runGeoSaveStepSTL(obj, name)
            % run gmsh via matlab
            pyPath = "C:\\Users\\Remote Server 3\\AppData\\Local\\Programs\\Python\\Python310\\python.exe";
            pyCodePath = "C:\\emdlab-win64\\py-files\\gmsh\\emdlab_gmsh_runGeoSaveStepSTL.py";

            [~,~] = system(char('"' + pyPath + '"' + " " + '"' + pyCodePath+ '"'));  

            stpPath = "C:\emdlab-win64\tmp\emdlab_gmsh_stpFile.step";
            stlPath = "C:\emdlab-win64\tmp\emdlab_gmsh_stlFile.stl";
            stpPathEntity = cd + "\" + name + ".step";
            stlPathEntity = cd + "\" + name + ".stl";
            copyfile(stpPath, stpPathEntity);
            copyfile(stlPath, stlPathEntity);

            obj.entities{end+1} = emdlab_g3d_entity;
            obj.entities{end}.name = name;
            obj.entities{end}.stpPath = stpPathEntity;
            obj.entities{end}.stlPath = stlPathEntity;
        end



        function readAndAddStepFile(obj)

            fid = fopen('emdlab_g3d_stepFile.step', 'r');
            str = fscanf(fid, '%s');

            obj.entities(end+1) = string(str);
            fclose(fid);

        end

        function subtractVolumes(obj, index1, index2)

            fid = fopen('emdlab_g3d_stepFile1.step', 'w');
            str = strsplit(obj.entities(1), ';');
            for i = 1:(numel(str)-1)
                fprintf(fid, '%s;\n', str{i});
            end
            fclose(fid);

            fid = fopen('emdlab_g3d_stepFile2.step', 'w');
            str = strsplit(obj.entities(2), ';');
            for i = 1:(numel(str)-1)
                fprintf(fid, '%s;\n', str{i});
            end
            fclose(fid);

            

            [~,~] = system('"C:\\Users\\Remote Server 3\\AppData\\Local\\Programs\\Python\\Python310\\python.exe" "C:\\emdlab-win64\\geometry\\step\\emdlab_stp_subtract.py"');
            
            obj.entities = [];
            obj.readAndAddStepFile
            
        end

       
        function addCylinde1r(obj, name, r)

            pyPath = "C:\\Users\\Remote Server 3\\AppData\\Local\\Programs\\Python\\Python310\\python.exe";
            stpPath = "C:\\emdlab-win64\\geometry\\step\\emdlab_g3d_stepFile.step";            
            pyCodePath = "C:\\emdlab-win64\\geometry\\emdlab_cq_script.py";

            fid = fopen("C:\emdlab-win64\geometry\emdlab_cq_script.py", 'w');
            fprintf(fid, 'import cadquery as cq\n');
            fprintf(fid, 'wp = cq.Workplane("XY")\n');
            fprintf(fid, "cylinder = wp.circle(%f).extrude(50)\n", r);
            fprintf(fid, 'cq.exporters.export(cylinder, "%s")\n', stpPath);        
            fclose(fid);
            
            [~,~] = system(char('"' + pyPath + '"' + " " + '"' + pyCodePath+ '"'));            
            copyfile(stpPath, cd + "\" + string(name) + ".step");
            

        end

        function openSteps(obj, varargin)

            fid = fopen("C:\emdlab-win64\tmp\emdlab_gmsh_geoFile.geo", 'w'); 
            fprintf(fid,'SetFactory("OpenCASCADE");\n');
            for i = 1:obj.Ne
                fprintf(fid,'Merge "%s";\n', obj.entities{i}.stpPath);
                fprintf(fid, 'Physical Volume(%d) = {%d};\n', i, i);
            end
            fprintf(fid,'Coherence;\n');
            fclose(fid);

        end

        function m = read_msh_file(obj)

            obj.openSteps;

            % run gmsh via matlab
            pyPath = "C:\\Users\\Remote Server 3\\AppData\\Local\\Programs\\Python\\Python310\\python.exe";
            pyCodePath = "C:\\emdlab-win64\\py-files\\gmsh\\emdlab_gmsh_runGeoSaveMsh3D.py";

%             [~,~] = system(char('"' + pyPath + '"' + " " + '"' + pyCodePath+ '"'));  

            % read generated mesh;
            emdlab_gmsh_mshFile;

            % get an instance of mesh data base
            m = emdlab_m3d_ttmdb;
            
            nodes = msh.POS(:,1:3);
            Np = size(nodes,1);

            % add faces
            for i = 1:numel(obj.entities)

                cl = msh.TETS(msh.TETS(:,5) == i, 1:4);
                index = unique(cl(:));
                index = sort(index);
                xpoints = nodes(index,:);
                pindex = zeros(Np,1);
                pindex(index) = 1:size(xpoints,1);
                cl = pindex(cl);

%                 p21 = xpoints(cl(:,2),:) - xpoints(cl(:,1),:);
%                 p31 = xpoints(cl(:,3),:) - xpoints(cl(:,1),:);
%                 index = (p21(:,1).*p31(:,2) - p21(:,2).*p31(:,1)) < 0;
%                 cl(index,:) = cl(index,[1,3,2]);

                m.addmz(obj.entities{i}.name, emdlab_m3d_ttmz(cl, xpoints));
%                 m.mzs.(obj.faces(i).tag).color = obj.faces(i).color;

            end

        end

        function subtract(obj, varargin)

            fid = fopen("C:\emdlab-win64\tmp\emdlab_gmsh_geoFile.geo", 'w'); 
            fprintf(fid,'SetFactory("OpenCASCADE");\n');            
            for i = 1:numel(varargin)
                fprintf(fid,'Merge "%s";\n', obj.getEntitySTPPathByName(varargin{i}));
            end
            fprintf(fid,'BooleanDifference{ Volume{1}; Delete; }{ Volume{2:%d}; Delete; }\n', numel(varargin));
            fprintf(fid, 'Coherence;\n');
            fclose(fid);

           % run gmsh via matlab
            pyPath = "C:\\Users\\Remote Server 3\\AppData\\Local\\Programs\\Python\\Python310\\python.exe";
            pyCodePath = "C:\\emdlab-win64\\py-files\\gmsh\\emdlab_gmsh_runGeoSaveStepSTL.py";

            [~,~] = system(char('"' + pyPath + '"' + " " + '"' + pyCodePath+ '"'));  

            stpPath = "C:\emdlab-win64\tmp\emdlab_gmsh_stpFile.step";
            stlPath = "C:\emdlab-win64\tmp\emdlab_gmsh_stlFile.stl";
            stpPathEntity = cd + "\" + varargin{1} + ".step";
            stlPathEntity = cd + "\" + varargin{1} + ".stl";
            copyfile(stpPath, stpPathEntity);
            copyfile(stlPath, stlPathEntity);

        end

       

        function showg(obj)

            system('python.exe C:\\emdlab-win64\\py-files\\vtk\\emdlab_vtk_showSTL.py')

        end

    end

end