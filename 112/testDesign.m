%DESIGN TO TEST AT SEA/STATIC 

clear 
clc
close all
warning off

%initial conditions and parameters
pr = [16:2:40]; %compressor pressure ratio
pr_burner = 0.95; %pressure recovery burner
etaPC = 0.90; %compressor polytropic efficiency

B = 0.3; %fan bypass ratio
etaF = 0.92; %adiabatic efficiency fan
etaFn = 0.99; %adiabatic efficiency fan nozzle
prf = 2; %fan pressure ratio

tmax = [1400:50:1800]; %t04 max combustion temp
M = 0.01; %flight mach
etaD = 0.95; %diffuser adiabatic efficiency
etaB = 0.97; %burner adiabatic efficiency
etaPT = 0.92; %turbine poltropic efficiency
etaN = 0.98; %nozzle adiabatic 
Q = 45000; %fuel heating
gamma = 1.4;
cp = 1.006;
R = 287;

%Ta = 220;
Ta = 288.15;
Pa = 101325;
rhoA = Pa/(R*Ta);

r = 0.9;
inletA = 3.1415*r^2;

PaSTD = 101325; %sea lvl std
TaSTD = 288.15; 

%preallocation 

%compressor inlet 
for i=1:length(tmax)
    for j=1:length(pr)
        gamma = 1.4;
        u = M*sqrt(gamma*R*Ta);
        rhoV = rhoA * u;
        mdotA = rhoV*inletA;
        isenRatio = (gamma-1)/2;
        T02 = Ta*(1+ isenRatio*M^2);
        P02 = Pa*(1+etaD*((T02/Ta)-1))^(gamma/(gamma-1));

        rhoVmax = 231.8*(P02/PaSTD)/sqrt(T02/TaSTD);
        %compressor outlet
        gamma = 1.37;
        P03 = P02*pr(j);
            %solve for adiabatic compressor efficency from polytropic
            syms etaAC
            fpc = (P03/P02)^(((gamma-1)/gamma)/etaPC) == 1+(1/etaAC)*((P03/P02)^((gamma-1)/gamma) - 1);
            etaC = vpasolve(fpc,etaAC);
        T03 = T02*(1+(1/etaC)*(pr(j)^((gamma-1)/gamma)-1));
        %burner and combustion
        gamma = 1.35;
        cp = ((R*gamma)/(gamma-1))/1000;
        f = ((tmax(i)/T03)-1)/(((Q*etaB)/(cp*T03))-(tmax(i)/T03));
        %turbine inlet 4
        P04 = P03*pr_burner;
        %fan outlet
        gamma = 1.4;
        P08 = P02*prf;
        T08 = T02*(1+(1/etaF)*(prf^((gamma-1)/gamma)-1));
        %fan nozzle exit
        gamma = 1.4;
        uef = sqrt(2*etaFn*(gamma/(gamma-1))*R*T08*(1-(Pa/P08)^((gamma-1)/gamma)));
        %turbine exit
        gamma = 1.33;
        T05 = tmax(i)-(T03-T02)-B*(T08-T02);
            %solve for adiabatic compressor efficency from polytropic
            syms etaAT
            fpt = (T05/tmax(i))^((gamma*etaPT)/(gamma-1)) == (1+((T05/tmax(i))-1)*etaAT)^(gamma/(gamma-1));
            etaT = vpasolve(fpt,etaAT);
        P05 = P04*(1-(1/etaT)*(1-(T05/tmax(i))))^(gamma/(gamma-1));
        %afterburner if applicable (nozzle inlet)change below
        P06 = P05;
        T06 = T05;
        %nozzle exit
        gamma = 1.36;
        ue = sqrt(2*etaN*(gamma/(gamma-1))*R*T06*(1-(Pa/P06)^((gamma-1)/gamma)));

        %thrust calculations
       rhoV = rhoA * u;
        mdotA = rhoV*inletA;
        STbare = ((1+f)*ue + B*uef - u*(1+B))/1000; %units of kNs/kg
        Tbare = STbare*mdotA;
        T = Tbare/(1.04+0.01*B^1.2); %thrust in KN
        ST = T/((mdotA)*(1+B));
        TSFC = f/ST;

%         TSFCbare = (f/((1+f)*ue + B*uef- u*(1+B)))*1000; %units of kg/kNs
%         STbare = TSFCbare/f;
%         Tbare = STbare*mdotA;
%         T = Tbare/(1.04+0.01*B^1.2);
%         ST = double(T/((mdotA)*(1+B)));
%         TSFC = double(f/ST);

%         ST = ((1+f)*ue - u)/1000; %units of kNs/kg
%         TSFC = f/((1+f)*ue-u)*1000; %units of kg/kNs

        arrayST(i,j) = double(ST);
        arrayTSFC(i,j) = double(TSFC);


        nth = (0.5*(1+f)*ue^2)-(0.5*u^2);
        nth = nth/(f*Q*1000); %1000 factor for kJ to J

        np = ST*u*1000; %multiply by 1000 to restore N
        np = np/((0.5*(1+f)*ue.^2)-(0.5*u^2));

        no = np*nth;
        hold on  
        plot(arrayST(1:end,j),arrayTSFC(1:end,j))
        
    end
    fprintf('Station %i/%i complete\n',i,length(tmax))    
    plot(arrayST(i,1:end),arrayTSFC(i,1:end))
end

xlabel('Specific Thrust kNs/kg')
ylabel('Thrust Specific Fuel Consumption kg/kNs')