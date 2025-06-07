
hold on
f = fopen('I:\My_Projetcs\flat_linear_pm_motor\Maxwell_sims\b.tab','r');

b = fscanf(f,'%f');

b = reshape(b,4,[])';


plot(b(:,1)*(180/pi)*(Nma/2)/(La/Nma),b(:,2))
plot(b(:,1)*(180/pi)*(Nma/2)/(La/Nma),b(:,3))
plot(b(:,1)*(180/pi)*(Nma/2)/(La/Nma),b(:,4))

% legend('MATLAB Program', 'ANSYS Maxwell')