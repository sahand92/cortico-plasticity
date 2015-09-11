[X, Y] = meshgrid((0:999) / 1000, (0:999) / 1000);
 p = [0.1, 0.2]
ineq1 = Y < p(2) * (1 - p(1));             %# First inequation
ineq2 = X < p(1) * (1 - (Y / (1 - p(1)))); %# Second inequation
both = ineq1 & ineq2;                      %# Intersection of both inequations

figure, hold on
c = 1:3;                                   %# Contour levels
contourf(c(1) * ineq1, [c(1), c(1)], 'b')

contourf(c(2) * ineq2, [c(2), c(2)], 'g')  %# Fill area for second inequation
contourf(c(3) * both, [c(3), c(3)], 'r')   %# Fill area for both inequations
legend('First', 'Second', 'Both')
alpha(0.5)