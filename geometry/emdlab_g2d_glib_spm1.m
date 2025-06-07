function emdlab_g2d_glib_spm1(g, Dsh, ORD, p, dm, embrace)

% input arguments check
arguments

    g (1,1) emdlab_g2d_db
    Dsh (1,1) double {mustBePositive}
    ORD (1,1) double {mustBePositive}
    p (1,1) double {mustBePositive,mustBeInteger}
    dm (1,1) double {mustBePositive}
    embrace (1,1) double {mustBePositive}

end

alpha_p = 2*pi/p;

o = emdlab_g2d_point(0,0);
p1 = emdlab_g2d_point(Dsh/2,0);
p2 = emdlab_g2d_point(ORD/2-dm,0);
p3 = emdlab_g2d_point(ORD/2,0);
p4 = p2.getRotateAroundOrigin(embrace*alpha_p/2);
p5 = p3.getRotateAroundOrigin(embrace*alpha_p/2);
p6 = p1.getRotateAroundOrigin(alpha_p/2);
p7 = p2.getRotateAroundOrigin(alpha_p/2);
p8 = p3.getRotateAroundOrigin(alpha_p/2);

o = g.addPoint(o);
p1 = g.addPoint(p1);
p2 = g.addPoint(p2);
p3 = g.addPoint(p3);
p4 = g.addPoint(p4);
p5 = g.addPoint(p5);
p6 = g.addPoint(p6);
p7 = g.addPoint(p7);
p8 = g.addPoint(p8);

e1 = g.addSegment(p1,p2);
e2 = g.addSegment(p2,p3);
e3 = g.addSegment(p4,p5);
e4 = g.addSegment(p6,p7);
e5 = g.addSegment(p7,p8);

e6 = g.addArc(o,p1,p6,1);
e7 = g.addArc(o,p2,p4,1);
e8 = g.addArc(o,p4,p7,1);
e9 = g.addArc(o,p3,p5,1);
e10 = g.addArc(o,p5,p8,1);

l1 = g.addLoop(e1,1,e7,1,e8,1,e4,0,e6,0);
l2 = g.addLoop(e2,1,e9,1,e3,0,e7,0);
l3 = g.addLoop(e3,1,e10,1,e5,0,e8,0);

g.addFace('rotor', l1);
g.addFace('magnet', l2);
g.addFace('rap', l3);


end