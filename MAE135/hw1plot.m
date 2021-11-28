v = 0:1:300;
rho = 1.2;
p = 100000;
pRatio = 1+(0.5*rho*v.^2)/p;
densityRatio = pRatio.^(1/1.4);
plot(v,densityRatio)
xlabel('Velocity (m/s)');
ylabel('Density Ratio \rho2/\rho1');
title('Density Ratio vs Velocity, Incomp. Stagnation')