function emdlab_g2d_lib_ipm1(g, Dsh, ORD, p)

% input arguments check
arguments

    g (1,1) emdlab_g2d_db
    Dsh (1,1) double {mustBePositive}
    ORD (1,1) double {mustBePositive}
    p (1,1) double {mustBePositive,mustBeInteger}

end

alpha_p = 2*pi/p;

o = emdlab_g2d_point(0,0);
p1 = emdlab_g2d_point(Dsh/2,0);
p2 = emdlab_g2d_point(ORD/2,0);
p3 = p2.getRotateAroundOrigin(alpha_p/2);
p4 = p1.getRotateAroundOrigin(alpha_p/2);


o = g.addPoint(o);
p1 = g.addPoint(p1);
p2 = g.addPoint(p2);
p3 = g.addPoint(p3);
p4 = g.addPoint(p4);

e1 = g.addSegment(p1,p2);
e2 = g.addArc(o,p2,p3,1);
e3 = g.addSegment(p3,p4);
e4 = g.addArc(o,p4,p1,0);

l1 = g.addLoop(e1,1,e2,1,e3,e1,e4,1);
g.addFace('rotor', l1);

end