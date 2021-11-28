clear
close
clc

h1=2.492/2;
%-0.235*atan(-x+1.5)+0.25
%0.2*atan(-5*(x-0.75))+0.32
%0.047*atan(9*x-11.5)+0.198
%-0.49*atan(-x+3.6)+0.658
syms x
y = piecewise(x<0.5,0.5,1>x>0.5,-0.141*atan(13*x-9.7)+0.322, 1.5>x>1, 0.047*atan(9*x-11.5)+0.198, 2.5>x>1.5,0.25, 3.5>x>2.5, 0.151*atan(5*x-15)+0.4296,3.51>x>3.5,1.22/2 );
fplot(y,'LineWidth',2)

hold on
xlim([0,3.5])
ylim([0,0.75])
xlabel('Length(m)')
ylabel('Height(m)')
title('2D Section View of Wind Tunnel - Half Height')
