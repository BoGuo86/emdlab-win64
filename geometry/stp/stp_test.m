

g = emdlab_stp_MODEL;

% g.addRectangle(0,0,1,1);
% g.addRectangle(3,3,1,1);
% g.addCircle(0,0,2);
% 
% g.addPolygon([0,5,8,3,2],[0,0,6,8,6])
g.addBox(0,0,0,5,5,5)
g.writeSTEP('ali1.step');

g = emdlab_stp_MODEL;
g.addBox(2,2,2,5,5,5)
% Write the STEP file
g.writeSTEP('ali2.step');


