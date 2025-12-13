% EMDLAB: Electrical Machines Design Laboratory
% two-dimensional magnetic-transient start point
% a star point is defined to force KCL on a number of voltage-fed stranded coils

classdef emdlab_solvers_mt2d_star < handle & matlab.mixin.SetGet

    properties

        % coil arm names: all coil arms are solid
        coilArms (1,:) string;

        % star connection index
        sci (1,1) double;

        % coil indicies
        ci (1,:) double;

        % voltage on star point
        voltage (1,:) = 0;

    end

end