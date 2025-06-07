% developer: https://ComProgExpert.com
% rotate points on x-y plane
% p is a matrix [Np x 2]
% rotAngle: is the amount of rotation [rad]
% xc, yc: center point of the rotation

function newP = emdlab_g2d_rotatePoints(p, rotAngle, xc, yc)

if nargin < 3
    xc = 0;
    yc = 0;
end

newP = zeros(size(p,1),2);

p(:,1) = p(:,1) - xc;
p(:,2) = p(:,2) - yc;

newP(:,1) = p(:,1) * cos(rotAngle) - p(:,2) * sin(rotAngle);
newP(:,2) = p(:,1) * sin(rotAngle) + p(:,2) * cos(rotAngle);

newP(:,1) = newP(:,1) + xc;
newP(:,2) = newP(:,2) + yc;

end