% rotor & magnet
% inner rotor surface-mounted permenent magnet motor

function emdlab_g2d_lib_rm_spm51(g, yShift, dm, dry, wp, embrace, name1, name2, name3)

% input arguments check
arguments

    g (1,1) emdlab_g2d_db
    yShift (1,1) double;
    dm (1,1) double {mustBePositive}
    dry (1,1) double {mustBePositive}
    wp (1,1) double {mustBePositive}
    embrace (1,1) double {mustBePositive}
    name1 (1,:) char;
    name2 (1,:) char;
    name3 (1,:) char;

end

s1 = g.addSegmentByCoordinates(0,yShift-dm-dry,wp/2,yShift-dm-dry);
s2 = g.extendSegmentBySegment(s1, pi/2, dry);
s3 = g.extendSegmentBySegment(s2, pi/2, (1-embrace)*wp/2);
s4 = g.extendSegmentBySegment(s3, 0, embrace*wp/2);
s5 = g.extendSegmentBySegment(s4, pi/2, dry);

s6 = g.extendSegmentBySegment(s3, -pi/2, dm);
s7 = g.extendSegmentBySegment(s6, pi/2, embrace*wp/2);
s8 = g.extendSegmentBySegment(s7, pi/2, dm);

s9 = g.extendSegmentBySegment(s2, 0, dm);
s10 = g.extendSegmentBySegment(s9, pi/2, (1-embrace)*wp/2);


l1 = g.addLoop(s1,s2,s3,s4,s5);
l2 = g.addLoop(s6,s7,s8,-s4);
l3 = g.addLoop(s9,s10,-s6,-s3);

g.addFace(name1, l1);
g.addFace(name2, l2);
g.addFace(name3, l3);

% set mesh zone colors
g.setFaceColor(name1,200,200,200)
g.setFaceColor(name2,28,255,28)
g.setFaceColor(name3,0,255,255)

end