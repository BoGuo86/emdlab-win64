clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));
g = emdlab_stp_MODEL;


fid = fopen('emdlab_g3d_stepFile.step', 'r');
while true

    if feof(fid)
        fclose(fid);
        break;
    end

    str = fgetl(fid);
    if numel(str)
        if strcmpi(str(1),'#')
            g.addEntityByText(str);
        end
    end

end

g.showModel
% rotate3d



