clear
clc
close

aS = 0.8418;
a4 = 1:0.5:6;

ratio = a4/aS;
mach = [ 0.60025348, 0.34909868, 0.25305104,  0.19955547, 0.16505258,  0.14084836,  0.12289576,  0.10902919,  0.09799204,  0.08899501,  0.08151607];

plot(mach,ratio, 'LineWidth',2)
xlabel('Mach Number')
ylabel('Area Ratio (A4/A*)')
ylim([0.5,7.5])
title('Exit Area and Exit Mach Number')
