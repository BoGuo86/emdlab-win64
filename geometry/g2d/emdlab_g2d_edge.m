% developer: https://ComProgExpert.com
% 2D edge class
% an edge can be:
% a segment
% an arc
% a spline

classdef emdlab_g2d_edge  < handle & matlab.mixin.Copyable & emdlab_g2d_constants

    properties

        % pointer to edge type
        ptr (1,1);

        % the object tag
        tag (1,:) char = '';

        % tag of loops that used this edge for their constrcution
        tags (1,:) string;

    end

end