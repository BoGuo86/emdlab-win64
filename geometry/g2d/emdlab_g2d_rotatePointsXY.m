% EMDLAB: Electrical Machines Design Laboratory
% rotate points on x-y plane
% rotAngle: is the amount of rotation [rad]
% xc, yc: center point of the rotation

function [xr, yr] = emdlab_g2d_rotatePointsXY(x, y, rotAngle, xc, yc)

if nargin < 4
    xc = 0;
    yc = 0;
end

x = x - xc;
y = y - yc;

xr = x * cos(rotAngle) - y * sin(rotAngle);
yr = x * sin(rotAngle) + y * cos(rotAngle);

xr = xr + xc;
yr = yr + yc;

end