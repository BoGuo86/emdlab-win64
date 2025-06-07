

clc
clear

m = emdlab_TMZPC([1,2,3;1,3,4]',[0,0;1,0;1,1;0,1]');

m.setMeshDataBase;

m.showm