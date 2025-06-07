function [f,v] = emdlab_g2d_getRectangleFVInside(x1,y1,x2,y2)

% facets
f = [1,2;
    2,3;
    3,4;
    4,1];

% vertices
v = [x1,y1;
    x2,y1;
    x2,y2;
    x1,y2];

end