
function p = emdlab_g2d_ftlfaptc(xc, yc, r, x0, y0)

p = zeros(2,2);
p(1,1) = (r^2*x0 - r^2*xc - 2*x0*xc^2 + x0^2*xc + xc*y0^2 + xc*yc^2 + xc^3 - r*y0*(- r^2 + x0^2 - 2*x0*xc + xc^2 + y0^2 - 2*y0*yc + yc^2)^(1/2) + r*yc*(- r^2 + x0^2 - 2*x0*xc + xc^2 + y0^2 - 2*y0*yc + yc^2)^(1/2) - 2*xc*y0*yc)/(x0^2 - 2*x0*xc + xc^2 + y0^2 - 2*y0*yc + yc^2);
p(2,1) =  (r^2*x0 - r^2*xc - 2*x0*xc^2 + x0^2*xc + xc*y0^2 + xc*yc^2 + xc^3 + r*y0*(- r^2 + x0^2 - 2*x0*xc + xc^2 + y0^2 - 2*y0*yc + yc^2)^(1/2) - r*yc*(- r^2 + x0^2 - 2*x0*xc + xc^2 + y0^2 - 2*y0*yc + yc^2)^(1/2) - 2*xc*y0*yc)/(x0^2 - 2*x0*xc + xc^2 + y0^2 - 2*y0*yc + yc^2);

p(1,2) =  (r^2*y0 - r^2*yc + x0^2*yc + xc^2*yc - 2*y0*yc^2 + y0^2*yc + yc^3 + r*x0*(- r^2 + x0^2 - 2*x0*xc + xc^2 + y0^2 - 2*y0*yc + yc^2)^(1/2) - r*xc*(- r^2 + x0^2 - 2*x0*xc + xc^2 + y0^2 - 2*y0*yc + yc^2)^(1/2) - 2*x0*xc*yc)/(x0^2 - 2*x0*xc + xc^2 + y0^2 - 2*y0*yc + yc^2);
p(2,2) =  (r^2*y0 - r^2*yc + x0^2*yc + xc^2*yc - 2*y0*yc^2 + y0^2*yc + yc^3 - r*x0*(- r^2 + x0^2 - 2*x0*xc + xc^2 + y0^2 - 2*y0*yc + yc^2)^(1/2) + r*xc*(- r^2 + x0^2 - 2*x0*xc + xc^2 + y0^2 - 2*y0*yc + yc^2)^(1/2) - 2*x0*xc*yc)/(x0^2 - 2*x0*xc + xc^2 + y0^2 - 2*y0*yc + yc^2);

[~, index] = sort(p(:,1));
p = p(index,:);

end