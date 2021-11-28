clear
clc
close

aS = 0.8418;
a4 = 1:0.5:6;

ratio = a4/aS;
mach = [ 0.60025348, 0.34909868, 0.25305104,  0.19955547, 0.16505258,  0.14084836,  0.12289576,  0.10902919,  0.09799204,  0.08899501,  0.08151607];
velocity = mach*sqrt(1.4*287*293.15);
pfan = 0.5.*velocity.^3.*a4.*1.2041./1000;
plot(pfan,ratio, 'LineWidth',2)
xlabel('Fan Power (kW)')
ylabel('Area Ratio (A4/A*)')

title('Exit Area and Fan Power')
