
function v = interpolate_vB(vB, B)

N = length(B);
v = zeros(N,1);

index = (B>=vB.breaks(1)) & (B<=vB.breaks(end));
v(index) = ppval(vB, B(index));

index = B>vB.breaks(end);
c = vB.coefs(end,:);
tmp = vB.breaks(end) - vB.breaks(end-1);
x_tmp = vB.breaks(end);
y_tmp = polyval(c, x_tmp);
slope = 3*c(1)*tmp^2+2*c(2)*tmp+c(3);
slope = (ppval(vB, vB.breaks(end)) - ppval(vB, vB.breaks(end-1)))/(vB.breaks(end)-vB.breaks(end-1));
v(index) = y_tmp + slope * (B(index) - x_tmp);


end