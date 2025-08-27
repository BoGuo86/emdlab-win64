clc;
clear;
close all;

g3 = emdlab_g3d_db;

g = emdlab_g2d_db;
g.addCircleLoop(0,0,100);
g.addFace('square',1);
g.extrudeAndSaveStep(0,20);

g3.readAndAddStepFile;

% g = emdlab_g2d_db;
% g.addCircleLoop(0,0,70);
% g.addFace('square',1);
% g.extrudeAndSaveStep(0,20);
% 
% g3.readAndAddStepFile;
% g3.subtractVolumes;

g = emdlab_g2d_db;
g.addRectangleLoop(0,-5,200,10);
g.addFace('square',1);
g.extrudeAndSaveStep(5,30);

g3.readAndAddStepFile;
g3.subtractVolumes;

% g.addCircleLoop(0,0,40);
% g.addCircleLoop(0,0,50);
% g.addFace('stator',2,1);
% g.setMeshMaxLength(5);
% g.showFaces
