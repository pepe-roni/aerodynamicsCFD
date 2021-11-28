clear
close
clc

h1=2.492/2;
%-1.48*atan(x)+2.41
%23.29*atan(x)-31.9

syms x
y = piecewise(x<1,h1,3.5>x>1,0.415*atan(-x+2.25)+0.873, 5.5>x>3.5, 0.5, 10.4>x>5.5,-0.825*atan(-x+7.95)+1.475, x>10.4, 4.901/2);
fplot(y,'LineWidth',2)

hold on
xlim([0,10.45])
ylim([0,2.5])
xlabel('Length(m)')
ylabel('Height(m)')
title('2D Section View of Wind Tunnel - Half Height')
