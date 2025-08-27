function [f,v] = emdlab_g2d_getRectangleFVOutside(x1,y1,x2,y2)

% facets
f = [1,4;
    4,3;
    3,2;
    2,1];

% vertices
v = [x1,y1;
    x2,y1;
    x2,y2;
    x1,y2];

end