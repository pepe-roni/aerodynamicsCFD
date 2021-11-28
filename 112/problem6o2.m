clear 
clc
mdot  = 100;
R =287;
Po = 16.186e3;
To = 257.964;
A = 3.07;
syms M P T
eq1 = P == Po/((1+((1.4-1)/2)*M^2)^(1.4/(1.4-1)));
eq2 = T == To/((1+((1.4-1)/2)*M^2));
eq3 = mdot/A == (P/(R*T))*M*sqrt(1.4*R*T);
sol = vpasolve([eq1,eq2, eq3],[M, P, T])
sol.M
sol.P
sol.T