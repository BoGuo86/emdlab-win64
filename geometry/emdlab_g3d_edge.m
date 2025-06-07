% developer: https://ComProgExpert.com
% 3D edge class
% an edge can be:
% a segment
% an arc
% a spline

classdef emdlab_g3d_edge  < handle & matlab.mixin.Copyable & emdlab_g2d_constants

    properties

        % pointer to edge type
        ptr (1,1);

        % the object tag
        tag (1,:) char = '';

    end

end