clear
clc
close

aS = 0.2058;
a4 = 0.4:0.2:1.6;

ratio = a4/aS;
mach = [ 0.3166, 0.20346608,   0.15091516,  0.12013145,   0.09984226,   0.08544259,   0.07468515];
velocity = mach*sqrt(1.4*287*324.55);
mdot = 1.088.*a4.*velocity;
pfan = mdot.*1006.1.*(327.15-293.15)/(10^-6);
plot(pfan,ratio, 'LineWidth',2)
xlabel('Fan Power (MW)')
ylabel('Area Ratio (A4/A*)')

title('Exit Area and Fan Power')
