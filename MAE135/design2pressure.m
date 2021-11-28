clear
close
clc

p01 = 101.325;
p04 = 104.191;
p02 = 144.509;
p03 = 144.509;

m1 = 0.24;
m2 = 0.172;
m3 = 2;
m3s = 0.5774;
m4 = 0.2;

a1 = -0.40971;
b1 = 3.95647;
c1 = -5.6846;
d1 = 3.64781;
f1 = -0.475119;

a2 = 0.0158038;
b2 = -9.49624;
c2 = 74.9795;
d2 = -112.349;
f2 = 47.85;

a3 = -4.714*10^-7;
b3 = 0.238367;
c3 = -10.3285;
d3 = 36.7026;
f3 = -35.9135;

syms x
y = piecewise(x<0.25,m1,0.25<x<0.5,m2,0.98>x>0.5, a1*x^11+b1*x^4+c1*x^2+d1*x+f1,1.18>x>0.98,x^2, 1.495>x>1.18, a2*x^12+b2*x^4+c2*x^2+d2*x+f2, 2.5>x>1.5, m3, x>2.5, a3*x^12+b3*x^4+c3*x^2+d3*x+f3);
%fplot(y,'LineWidth',2)
pRatio = (1+((1.4-1)/2)*y^2)^(1.4/(1.4-1))
p1 = p01/pRatio;
p2 = p02/pRatio;
pstar = p02/pRatio;
p3 = p03/pRatio;
p4 = p04/pRatio;
p = piecewise(x<0.25,p1,0.25<x<0.5,p2,1>x>0.5,pstar, 1.5>x>1.0,pstar, 2.5>x>1.5,p3, x>2.5,p4);
fplot(p)

hold on
xlim([0,3.5])
ylim([0,145])
ylabel('Pressure (kPa)')
xlabel('Length(m)')
title('Pressure Through the Length of the Tunnel')