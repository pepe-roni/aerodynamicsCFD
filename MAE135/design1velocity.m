clear
close
clc

p01 = 101.325;
p04 = 102.034;
p02 = 102.034;
p03 = 102.034;

t01 = 293.15;
t02 = 293.735;
t03 = 293.735;
t04 = 293.735;

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
y = piecewise(x<0.5,m1,0.5<x<1,m2,3.5>x>1, a1*x^3+b1*x^2+c1*x+d1, 5.51>x>3.5, 0.6, 10.4>x>5.51, a2*x^9+b2*x^5+c2*x^2+d2*x+g, x>10.4, m4);
%fplot(y,'LineWidth',2)
pRatio = (1+((1.4-1)/2)*y^2)^(1.4/(1.4-1))
p1 = p01/pRatio;
p2 = p02/pRatio;
p3 = p03/pRatio;
p4 = p04/pRatio;
p = piecewise(x<0.5,p1,0.5<x<1,p2,3.5>x>1, p3, 5.51>x>3.5, p3, 10.4>x>5.51, p4, x>10.4,p4);
tRatio = (1+((1.4-1)/2)*y^2)
t1 = t01/tRatio;
t2 = t02/tRatio;
t3 = t03/tRatio;
t4 = t04/tRatio;
t = piecewise(x<0.5,t1,0.5<x<1,t2,3.5>x>1, t3, 5.51>x>3.5, t3, 10.4>x>5.51, t4, x>10.4,t4);
rho = p/(287*t);
v = y * sqrt(1.4*287*t)
fplot(v)
hold on
xlim([0,10.3])
ylim([10,210]);

ylabel('Velocity(m/s)')
xlabel('Length(m)')
title('Velocity Through the Length of the Tunnel')