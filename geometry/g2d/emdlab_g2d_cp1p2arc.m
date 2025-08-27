function z = emdlab_g2d_cp1p2arc(x1, y1, x2, y2, x3, y3, varargin)

z = emdlab_g2d_arc(emdlab_g2d_point(x1, y1), emdlab_g2d_point(x2,y2), emdlab_g2d_point(x3,y3), varargin{:});

end