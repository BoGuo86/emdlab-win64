
clc

N = 3;
R = 175/2;

while true
  chk = (N/pi)*cos(pi/N)*sin(pi/N);
  if chk > 0.99
    break
  end
  N = N + 1;
end

disp(2*pi*R/N)


