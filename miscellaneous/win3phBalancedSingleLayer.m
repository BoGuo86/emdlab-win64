
function win3phBalancedSingleLayer(p,Q)

slot = (1:Q)';

Tu = (slot-1)*(p/2)*(360/Q);
Tu = rem(Tu,360);

disp(Tu)


end