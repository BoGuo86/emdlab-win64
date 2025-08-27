classdef emdlab_g2d_location
    enumeration
        % location of a point with respect to a segment
        ps_right
        ps_left
        ps_beyond
        ps_behind
        ps_between
        ps_origin
        ps_destination

        % location of a point with respect to a line
        pl_right
        pl_left
        pl_on

        % location of a point with respect to a circle
        pc_in
        pc_on
        pc_out
    end
end
