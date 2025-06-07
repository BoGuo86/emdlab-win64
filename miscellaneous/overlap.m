
hold on
f = fopen('I:\My_Projetcs\flat_linear_pm_motor\Maxwell_sims\By.tab','r');

b = fscanf(f,'%f');

b = reshape(b,2,[])';

plot(b(:,1),b(:,2))

legend('MATLAB Program', 'ANSYS Maxwell')