clear 
clc
close all
warning off

%initial conditions and parameters
rd = 0.7;
rb = 0.95;
rn = 0.98;
q = 45000000;
Ta = 200;
cp = 1004;
nb = 1;
gamma = 1.4;

tmax = [1640, 2500, 3000];

%design process

Ma = (1:0.1:8);
fRight = rd*rb*rn*(1+((gamma-1)/2)*Ma.^2).^(gamma/(gamma-1));
fRight = fRight.^((gamma-1)/gamma);
fRight = (fRight - 1)/((gamma-1)/2);
Me = sqrt(fRight);

T02 = Ta*(1+((gamma-1)/2)*Ma.^2);

f = zeros(1,length(Ma));
for i=1:length(tmax)
    colorstring = 'mbgrkc';
    for j=1:length(Ma)
        f(j) = ((tmax(i)/T02(j))-1)/(((q*nb)/(cp*T02(j)))-(tmax(i)/T02(j)));
    end
    Te = tmax(i)./(1+((gamma-1)/2)*Me.^2);
    ue = Me.*sqrt(gamma*287.*Te);
    u = Ma.*sqrt(gamma*287.*Ta);
    ST = (1+f).*ue - u;
    ST = ST./1000; %convert from kNs/g to kNs/kg
    TSFC = f./ST;
    
    
    hold on
    figure(1)
    yyaxis left
    title('Ramjet Thrust and Fuel Consumption');
    xlabel('Mach Number')
    ylabel('Specific Thrust, kN*s/kg')
    plot(Ma,ST,'-','Color', colorstring(i))
    ylim([0,1.200])
    yyaxis right
    ylabel('Thrust Specific Fuel Consumption,kg/kN*s')      
    plot(Ma(1:15+ST*120),TSFC(1:15+ST*120),'--','Color', colorstring(i+3)) %plot with varied TSFC
    ylim([0,0.200])
    legend('2000K - ST', '2500K - ST', '3000K - ST','2000K - TSFC', '2500K - TSFC', '3000K - TSFC')
  
    if(i==2)
        figure(2)
        hold on
        title('Ramjet Thrust and Efficiencies @Tmax = 2500K')
        xlabel('Mach Number')
        ylabel('Specific Thrust, kN*s/kg')
        plot(Ma,ST,'-')
        ylim([0,1.00])
        nth = (0.5*(1+f).*ue.^2)-(0.5*u.^2);
        nth = nth./(f*q);
        plot(Ma(1:45),nth(1:45))
        np = ST.*u*1000; %multiply by 1000 to restore N
        np = np./((0.5*(1+f).*ue.^2)-(0.5*u.^2));
        plot(Ma(1:45),np(1:45))
        no = np.*nth;
        plot(Ma(1:45),no(1:45));
        legend('2500K - ST','Propulsive Efficiency','Thermal Efficiency', 'Overall Efficiency');
    end 
end



