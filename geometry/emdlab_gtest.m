% initialization
clc;
clear;
close all;
addpath(genpath('C:\emdlab-win64'));
% define geometry data base
g = emdlab_g2d_db;

ri = 24/2;
ro = 74/2;
p = 4;
dm = 3;
wrrib = 1.5;
wtrib = 1;
alpha_p = 2*pi/p;

l = emdlab_g2d_line(emdlab_g2d_point(0,0),emdlab_g2d_point(cos(alpha_p/2),sin(alpha_p/2)));

s = g.addSegmentByCoordinates(0,wrrib/2,ri,wrrib/2);
s = g.extentSegmentBySegment(s,1,pi/3,10);
s = g.extentSegmentByTangentArc(s,1,dm/2,pi);
% s = g.extentArcByTangentSegment(s,1,g.edges(7).ptr.getLength);
% g.points(6).x = g.points(13).x;

g.addCircleLoop(0,0,ri);
g.addCircleLoop(0,0,ro);

g.setMeshMaxLength(1);
g.showSketch

% g.addSegmentByCoordinates(ri,0,ro,0);
% g.addArcByCoordinatesCPA(0,0,ro,0,alpha_p/2);
% g.addArcByCoordinatesCPA(0,0,ri,0,alpha_p/2);
% g.addSegment(4,5);
% 

% 
% s = g.extentSegmentByTangentArc(5,0,getDistanceFromOrigin(g.points(13)-g.points(6))/2,pi);
% 
% g.addSegment(10,12);
% g.addSegment(13,9);
% % g.extentSegmentBySegment(5,1,pi/3,16);
% % g.extentSegmentBySegment(5,0,pi/3,16);
% % g.addSegment(8,9);
% 
% g.setMeshMaxLength(0.2);
% g.showSketch(0);
% 
% l1 = g.addLoop(1,1,2,1,4,1,3,0);
% 
% l2 = g.addLoop(5,1,6,1,7,1,8,1,9,1,10,1);
% 
% l3 = g.addLoop(11,1,9,1,12,1,7,1);
% 
% l4 = g.addLoop(11,0,8,1);
% l5 = g.addLoop(10,1,5,1,6,1,12,0);
% 
% g.addFace('f',l1,l2);
% g.addFace('m',l3);
% 
% g.addFace('a1',l4);
% g.addFace('a2',l5);
% 
% 
% 
% e = g.getEdgeHandleByIndex(10);
% e.ptr.setMaxLength(2);
% % g.showSketch(0);
% % 
% 
% 
% m.setMeshZoneColor('f',200,200,200)
% m.setMeshZoneColor('m',28,255,28)
% m.joinMeshZones('a', 'a1', 'a2')
% m.setMeshZoneColor('a',0,255,255)
% 
% m.aux_cmrjmzx('f',p,1)
% m.aux_cmrmzx('m',p,1)
% m.aux_cmrjmzx('a',p,1)
% m.showmzs

% g.addSegmentByCoordinates(ri,0,ro,0);
% g.addArcByCoordinatesCPA(0,0,ro,0,alpha_p/2);
% g.addArcByCoordinatesCPA(0,0,ri,0,alpha_p/2);
% g.addSegment(4,5);
% 
% s = g.addSegmentByCoordinates(0,0.5,35,0.5);
% s = g.extentSegmentByTangentArc(s,1,1,1*pi/3);
% s = g.extentArcByTangentSegment(s,1,13);
% s = g.extentSegmentByTangentArc(s,1,2,pi);
% s = g.extentArcByTangentSegment(s,1,g.edges(7).ptr.getLength);
% g.points(6).x = g.points(13).x;
% 
% s = g.extentSegmentByTangentArc(5,0,getDistanceFromOrigin(g.points(13)-g.points(6))/2,pi);
% 
% g.addSegment(10,12);
% g.addSegment(13,9);
% % g.extentSegmentBySegment(5,1,pi/3,16);
% % g.extentSegmentBySegment(5,0,pi/3,16);
% % g.addSegment(8,9);
% 
% g.setMeshMaxLength(0.2);
% g.showSketch(0);
% 
% l1 = g.addLoop(1,1,2,1,4,1,3,0);
% 
% l2 = g.addLoop(5,1,6,1,7,1,8,1,9,1,10,1);
% 
% l3 = g.addLoop(11,1,9,1,12,1,7,1);
% 
% l4 = g.addLoop(11,0,8,1);
% l5 = g.addLoop(10,1,5,1,6,1,12,0);
% 
% g.addFace('f',l1,l2);
% g.addFace('m',l3);
% 
% g.addFace('a1',l4);
% g.addFace('a2',l5);
% 
% 
% 
% e = g.getEdgeHandleByIndex(10);
% e.ptr.setMaxLength(2);
% % g.showSketch(0);
% % 

% 
% m.setMeshZoneColor('f',200,200,200)
% m.setMeshZoneColor('m',28,255,28)
% m.joinMeshZones('a', 'a1', 'a2')
% m.setMeshZoneColor('a',0,255,255)
% 
% m.aux_cmrjmzx('f',p,1)
% m.aux_cmrmzx('m',p,1)
% m.aux_cmrjmzx('a',p,1)
% m.showmzs

