clear
close
clc

m1 = 0.197;
m2 = 0.2;
m3 = 0.6;
m4 = 0.1;

a1 = 0.00369453;
b1 = 0.0420651;
c1 = -0.0919563;
d1 = 0.248209;

a2 = -5.0824*10^-10;
b2 = 0.0000187684;
c2 = -0.0193177;
d2 = -0.0947946;
g = 1.61616;

syms x
y = piecewise(x<0.5,m1,0.5<x<1,m2,3.5>x>1, a1*x^3+b1*x^2+c1*x+d1, 5.5>x>3.5, 0.6, 10.4>x>5.5, a2*x^9+b2*x^5+c2*x^2+d2*x+g, x>10.4, m4);
fplot(y,'LineWidth',2)

hold on
xlim([0,10.3])
ylim([0,0.8])
ylabel('Mach Number')
xlabel('Length(m)')
title('Mach Number Through the Length of the Tunnel')
