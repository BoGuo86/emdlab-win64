function z = emdlab_g2d_cparc(xc, yc, x1, y1, a, varargin)

c = emdlab_g2d_point(xc, yc);
p0 = emdlab_g2d_point(x1, y1);
z = emdlab_g2d_arc(c, p0, p0.getRotateAroundPoint(c, a), varargin{:});

end