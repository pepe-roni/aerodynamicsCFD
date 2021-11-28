clear
clc

pRatio = [1.001:0.25:20];

syms M ptRatio Ms

for i=1:length(pRatio)
    eq1 = pRatio(i) == ((2*1.4*M^2)-(1.4-1))/(1.4+1);
    Mcalc = vpasolve(eq1,M); %solve isen eq for M
    M0(i) = double(Mcalc(2)); %take only the positive values, load array
    
    eq2 = Ms^2 == ((1.4-1)*M0(i)^2 + 2)/((2*1.4*M0(i)^2)-(1.4-1)); %shock eq
    Mshock = vpasolve(eq2,Ms);
    MS(i) = double(Mshock(2));
    
    eq3 = ptRatio ==  (((1.4+1)*M0(i)^2)/((1.4-1)*M0(i)^2 +2))^(1.4/(1.4-1))*((1.4+1)/(2*1.4*M0(i)^2-(1.4-1)))^(1/(1.4-1));
    ptCalc = vpasolve(eq3, ptRatio);
    PT(i) = double(ptCalc);
end

figure(1)
plot(pRatio,MS)
title('Mach Shock Number as a function of Static Pressure Ratio')
xlabel('Static Pressure Ratio, P2/Pa')
ylabel('Mach Shock Number')

figure(2)
plot(pRatio,PT)
title('Total Pressure Ratio as a function of Static Pressure Ratio')
xlabel('Static Pressure Ratio, P2/Pa')
ylabel('Total Pressure Ratio, Po2/Poa')