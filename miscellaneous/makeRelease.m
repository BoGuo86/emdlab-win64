
if exist('EMDLAB_Release','file')
    rmdir('EMDLAB_Release','s')
end

mkdir('EMDLAB_Release')
mkdir('EMDLAB_Release\g2d')
mkdir('EMDLAB_Release\mex')
mkdir('EMDLAB_Release\pcode')
mkdir('EMDLAB_Release\mex\glib')
mkdir('EMDLAB_Release\mex\ext_functions')
mkdir('EMDLAB_Release\mec')
mkdir('EMDLAB_Release\pec')

mFileDir = {'Important_Functions','Mesh_Classes','Important_Classes','g2d'};

for i = 1:numel(mFileDir)
    myFiles = what([cd,'\',mFileDir{i}]);
    for j = 1:numel(myFiles.m)
        fname = myFiles.m{j};
        pcode(fname);
        if strcmp(fname(1:5),'mlib_')
            movefile([fname(1:end-2),'.p'],'EMDLAB_Release\mlib');
        elseif strcmp(fname(1:4),'g2d_')
            movefile([fname(1:end-2),'.p'],'EMDLAB_Release\g2d');
        else
            movefile([fname(1:end-2),'.p'],'EMDLAB_Release\pcode');
        end
    end
end

mFileDir = {'MEX'};

for i = 1:numel(mFileDir)
    myFiles = what([cd,'\',mFileDir{i}]);
    for j = 1:numel(myFiles.mex)
        if strcmp(myFiles.mex{j}(1:5),'glib_')
            copyfile(['MEX\',myFiles.mex{j}],'EMDLAB_Release\mex\glib');
        elseif strcmp(myFiles.mex{j}(1:4),'ext_')
            copyfile(['MEX\',myFiles.mex{j}],'EMDLAB_Release\mex\ext_functions');
        else
            copyfile(['MEX\',myFiles.mex{j}],'EMDLAB_Release\mex');
        end
    end
end


list = dir('Solvers_Classes');
list = list([list.isdir]);
list = list(~ismember({list.name},{'.' '..'}));

mkdir('EMDLAB_Release\solvers_classes')

for i = 1:numel(list)
    
    tmp = ['EMDLAB_Release\solvers_classes\',list(i).name];
    mkdir(tmp);
    
    myFiles = what([cd,'\solvers_classes\',list(i).name]);
    for j = 1:numel(myFiles.m)
        fname = myFiles.m{j};
        pcode(fname);
        movefile([fname(1:end-2),'.p'],tmp);
    end
end

list = dir('mlib');
list = list([list.isdir]);
list = list(~ismember({list.name},{'.' '..'}));

mkdir('EMDLab_Release\mlib')

for i = 1:numel(list)
    
    tmp = ['EMDLAB_Release\mlib\',list(i).name];
    mkdir(tmp);
    
    myFiles = what([cd,'\mlib\',list(i).name]);
    for j = 1:numel(myFiles.m)
        fname = myFiles.m{j};
        pcode(fname);
        movefile([fname(1:end-2),'.p'],tmp);
    end
end

list = dir('Examples');
list = list([list.isdir]);
list = list(~ismember({list.name},{'.' '..'}));

mkdir('EMDLAB_Release\examples')

for i = 1:numel(list)
    
    tmp = ['EMDLAB_Release\examples\',list(i).name];
    mkdir(tmp);
    
    myFiles = what([cd,'\examples\',list(i).name]);
    for j = 1:numel(myFiles.m)
        fname = myFiles.m{j};
        copyfile(['examples\',list(i).name,'\',fname],tmp);
    end
end

mkdir('EMDLAB_Release\materials')
tmp = dir('materialsdata');

for i = 1:numel(tmp)
    str = tmp(i).name;
    if length(str)>4
    if strcmpi(str(end-3:end), '.dat') || strcmpi(str(end-1:end), '.m')
    copyfile(['materialsdata\', tmp(i).name],'emdlab_release\materials')
    end
    end
end
