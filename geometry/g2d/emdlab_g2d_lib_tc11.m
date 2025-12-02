% stator tooth & coil
% rectangular slot & coil

function emdlab_g2d_lib_tc11(g, ID, OD, Ns, ws, ds, d0, name1, name2, name3)

% input arguments check
arguments

    g (1,1) emdlab_g2d_db
    ID (1,1) double {mustBePositive}
    OD (1,1) double {mustBePositive}
    Ns (1,1) double {mustBePositive,mustBeInteger}
    ws (1,1) double {mustBePositive}
    ds (1,1) double {mustBePositive}
    d0 (1,1) double {mustBePositive}
    name1 (1,:) char
    name2 (1,:) char
    name3 (1,:) char

end

gamma_so = 2*asin(ws/ID);
alpha_s = 2*pi/Ns;

o = emdlab_g2d_point(0,0);
p1 = emdlab_g2d_point((ID/2) * cos(gamma_so/2), (ID/2) * sin(gamma_so/2));
p2 = emdlab_g2d_point(ID/2 + ds, ws/2);
p3 = emdlab_g2d_point(ID/2 + ds, 0);
p4 = emdlab_g2d_point(OD/2, 0);
p5 = p4.getRotateAroundPoint(o,alpha_s/2);
p6 = p5.getMoveRadial(-OD/2+ID/2);
p7 = emdlab_g2d_point(ID/2, 0);
p8 = p7.getMoveX(d0);
p9 = p8.getMoveX(ds-2*d0);
p10 = p9.getMoveY(ws/2-d0);
p11 = p8.getMoveY(ws/2-d0);

o = g.addPoint(0,0);
p1 = g.addPoint(p1);
p2 = g.addPoint(p2);
p3 = g.addPoint(p3);
p4 = g.addPoint(p4);
p5 = g.addPoint(p5);
p6 = g.addPoint(p6);
p7 = g.addPoint(p7);
p8 = g.addPoint(p8);
p9 = g.addPoint(p9);
p10 = g.addPoint(p10);
p11 = g.addPoint(p11);

e1 = g.addSegment(p1,p2);
e2 = g.addSegment(p2,p3);
e3 = g.addSegment(p3,p4);
e4 = g.addArc(o,p4,p5,1);
e5 = g.addSegment(p5,p6);
e6 = g.addArc(o,p6,p1,0);e6 = g.addArc(o,p6,p1,0);

l1 = g.addLoop(e1,e2,e3,e4,e5,e6);

e7 = g.addSegment(p8,p9);
e8 = g.addSegment(p9,p10);
e9 = g.addSegment(p10,p11);
e10 = g.addSegment(p11,p8);

e11 = g.addSegment(p7,p8);
e12 = g.addSegment(p9,p3);
e13 = g.addArc(o,p1,p7,0);

l2 = g.addLoop(e7,e8,e9,e10);
l3 = g.addLoop(e11,-e10,-e9,-e8,e12,-e2,-e1,e13);

g.addFace(name1, l1);
g.addFace(name2, l2);

g.addFace(name3, l3);

g.setFaceColor(name1,200,200,200);
g.setFaceColor(name2,255,137,39);
g.setFaceColor(name3,0,255,255);

end