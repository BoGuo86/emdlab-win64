
a1 = G2DSEGMENT([0,0],[1,0]);
a2 = G2DARC([0,0],[1,0],[0,1]);
a = mexTest;
a = G2DARC(a(1,:),a(2,:),a(3,:));
hold all
a1.show('a1');
a2.show('a2');
a.show('a');