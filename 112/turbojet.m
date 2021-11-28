clear 
clc
close all
warning off

%initial conditions and parameters
pr = 30; %compressor pressure ratio
pr_burner = 1;

tmax = 1700;
M = 0.85;
etaD = 0.97;
etaC = 0.85;
etaB = 1;
etaT = 0.9;
etaN = 0.98;
Q = 45000;
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
%turbine inlet
P04 = P03(pr_burner);
%turbine exit
T05 = tmax - (T03-T02);
P05 = P04*(1-(1/etaT)*(1-(T05/tmax)))^(gamma/(gamma-1));
%afterburner if applicable (nozzle inlet)change below
P06 = P05;
T06 = T05;
%nozzle exit
ue = sqrt(2*etaN*(gamma/(gamma-1))*R*T06*(1-(Pa/P06)^((gamma-1)/gamma)));

%thrust calculations
ST = ((1+f)*ue - u)/1000 %units of kNs/kg
TSFC = f/((1+f)*ue-u)*1000 %units of kg/kNs

nth = (0.5*(1+f)*ue^2)-(0.5*u^2);
nth = nth/(f*Q*1000) %1000 factor for kJ to J

np = ST*u*1000; %multiply by 1000 to restore N
np = np/((0.5*(1+f)*ue.^2)-(0.5*u^2))

no = np*nth