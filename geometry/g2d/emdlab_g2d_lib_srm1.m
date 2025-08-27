function emdlab_g2d_lib_srm1(g, gv_Ns, gv_Nr, gv_Dsh, gv_ISD, gv_OSD, gv_gap, gv_beta_r, gv_beta_s, gv_wry, gv_wsy)

% dependent variables
gv_wrt = 2 * (gv_ISD/2-gv_gap) * sin(gv_beta_r/2);
gv_gammar = asin(gv_wrt*0.5/(gv_Dsh/2+gv_wry));
gv_alpharp = 2*pi/gv_Nr;
gv_wst = 2 * (gv_ISD/2) * sin(gv_beta_s/2);
gv_gammas = asin(gv_wst*0.5/(gv_OSD/2-gv_wsy));
gv_alphasp = 2*pi/gv_Ns;

% adding points
p1t = g.addPoint(0,0);
p2t = g.addPoint(gv_Dsh/2,0);
p3t = g.addPoint(gv_ISD/2 - gv_gap,0);
p4t = g.addPoint((gv_ISD/2 - gv_gap)*cos(gv_beta_r/2),(gv_ISD/2 - gv_gap)*sin(gv_beta_r/2));
p5t = g.addPoint((gv_Dsh/2+gv_wry)*cos(gv_gammar),(gv_Dsh/2+gv_wry)*sin(gv_gammar));
p6t = g.addPoint((gv_Dsh/2+gv_wry)*cos(gv_alpharp/2),(gv_Dsh/2+gv_wry)*sin(gv_alpharp/2));
p7t = g.addPoint((gv_Dsh/2)*cos(gv_alpharp/2),(gv_Dsh/2)*sin(gv_alpharp/2));
p8t = g.addPoint((gv_ISD/2 - gv_gap)*cos(gv_alpharp/2),(gv_ISD/2 - gv_gap)*sin(gv_alpharp/2));
p9t = g.addPoint(gv_ISD/2,0);
p10t = g.addPoint(gv_OSD/2,0);
p11t = g.addPoint((gv_OSD/2)*cos(gv_alphasp/2),(gv_OSD/2)*sin(gv_alphasp/2));
p12t = g.addPoint((gv_OSD/2-gv_wsy)*cos(gv_alphasp/2),(gv_OSD/2-gv_wsy)*sin(gv_alphasp/2));
p13t = g.addPoint((gv_OSD/2-gv_wsy)*cos(gv_gammas),(gv_OSD/2-gv_wsy)*sin(gv_gammas));
p14t = g.addPoint((gv_ISD/2)*cos(gv_beta_s/2),(gv_ISD/2)*sin(gv_beta_s/2));
p15t = g.addPoint((gv_ISD/2)*cos(gv_alphasp/2),(gv_ISD/2)*sin(gv_alphasp/2));

% adding edges
e1t = g.addSegment(p2t, p3t);
e2t = g.addArc(p1t, p3t, p4t, 1);
e3t = g.addSegment(p4t, p5t);
e4t = g.addArc(p1t, p5t, p6t, 1);
e5t = g.addSegment(p6t, p7t);
e6t = g.addArc(p1t, p7t, p2t, 0);
e7t = g.addArc(p1t, p4t, p8t, 1);
e8t = g.addSegment(p8t, p6t);
e9t = g.addSegment(p9t, p10t);
e10t = g.addArc(p1t, p10t, p11t, 1);
e11t = g.addSegment(p11t, p12t);
e12t = g.addArc(p1t, p12t, p13t, 0);
e13t = g.addSegment(p13t, p14t);
e14t = g.addArc(p1t, p14t, p9t, 0);
e15t = g.addArc(p1t, p15t', p14t, 0);
e16t = g.addSegment(p12t, p15t);

% adding loops
l1t = g.addLoop(e1t,e2t,e3t,e4t,e5t,e6t);
l2t = g.addLoop(-e3t,e7t,e8t,-e4t);
l3t = g.addLoop(e9t,e10t,e11t,e12t,e13t,e14t);
l4t = g.addLoop(-e13t,-e12t,e16t,e15t);

% adding faces
g.addFace('Rotor', l1t);
g.addFace('RotorAP', l2t);
g.addFace('Stator', l3t);
g.addFace('sca', l4t);

end