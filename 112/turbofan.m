clear 
clc
close all
warning off

%initial conditions and parameters
pr = 30; %compressor pressure ratio
pr_burner = 1; %pressure recovery burner
etaC = 0.85; %compresor adiabatic efficiency

B = 8 ; %fan bypass ratio
etaF = 0.85; %adiabatic efficiency fan
etaFn = 0.97; %adiabatic efficiency fan nozzle
prf = 1.5; %fan pressure ratio

tmax = 1800; %t04 max combustion temp
M = 0.85; %flight mach
etaD = 0.97; %diffuser adiabatic efficiency
etaB = 1; %burner adiabatic efficiency
etaT = 0.9; %trubine adiabatic efficiency
etaN = 0.98; %nozzle adiabatic 
Q = 45000; %fuel heating
gamma = 1.4;
cp = 1.006;
R = 287;

%Ta = 220;
Ta = 216.7;
Pa = 18.75e3;

%compressor inlet
u = M*sqrt(gamma*R*Ta);
isenRatio = (gamma-1)/2;
T02 = Ta*(1+ isenRatio*M^2);
P02 = Pa*(1+etaD*((T02/Ta)-1))^(gamma/(gamma-1));
%compressor outlet
P03 = P02*pr;
T03 = T02*(1+(1/etaC)*(pr^((gamma-1)/gamma)-1));
%burner and combustion
f = ((tmax/T03)-1)/(((Q*etaB)/(cp*T03))-(tmax/T03));
%turbine inlet 4
P04 = P03(pr_burner);
%fan outlet
P08 = P02*prf;
T08 = T02*(1+(1/etaF)*(prf^((gamma-1)/gamma)-1));
%fan nozzle exit
uef = sqrt(2*etaFn*(gamma/(gamma-1))*R*T08*(1-(Pa/P08)^((gamma-1)/gamma)));
%turbine exit
T05 = tmax-(T03-T02)-B*(T08-T02);
P05 = P04*(1-(1/etaT)*(1-(T05/tmax)))^(gamma/(gamma-1));
%afterburner if applicable (nozzle inlet)change below
P06 = P05;
T06 = T05;
%nozzle exit
ue = sqrt(2*etaN*(gamma/(gamma-1))*R*T06*(1-(Pa/P06)^((gamma-1)/gamma)));

%thrust calculations
ST = ((1+f)*ue + B*uef - u*(1+B))/1000 %units of kNs/kg
TSFC = (f/((1+f)*ue + B*uef- u*(1+B)))*1000 %units of kg/kNs

nth = (0.5*(1+f)*ue^2)-(0.5*u^2);
nth = nth/(f*Q*1000) %1000 factor for kJ to J

np = ST*u*1000; %multiply by 1000 to restore N
np = np/((0.5*(1+f)*ue.^2)-(0.5*u^2))

no = np*nth