clear
clc
close

aS = 0.2058;
a4 = 0.4:0.2:1.6;

ratio = a4/aS;
mach = [ 0.3166, 0.20346608,   0.15091516,  0.12013145,   0.09984226,   0.08544259,   0.07468515];

plot(mach,ratio, 'LineWidth',2)
xlabel('Mach Number')
ylabel('Area Ratio (A4/A*)')
ylim([0.5,8])
title('Exit Area and Exit Mach Number')
