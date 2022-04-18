clc
clear

%specs
range = 4000; %range in miles
Vao = 135; %landing approach speed in knots
Mcruise = 0.82; %cruise mach number
airfoilType = 'c'; 
sweep = 35;
AR = 8;
sigma = 0.953; %density ratio
PR = 0.2360; %pressure ratio

%insert selected airfoil t/c curve
airfoilSweepCurve = @(Md) 0.0071*Md^2 - 0.4786*Md + 0.4898; %conventional 35deg

%clmax estimates
CLestimateLand = @(x) 109.22*x^3 - 67.377*x^2 + 16.45*x + 2.002; %landing Cl from "c"
CLestimateTO = @(x) 84.325414*x^3 - 65.343054*x^2 + 17.178855*x + 1.029074;  %TO Cl from "c"

err = 1; Cl = 0.55; %inital guess

    if airfoilType == 'c'
        dMdiv = -0.2044*Cl^2 - 0.1122*Cl + 0.1236;
    elseif airfoilType == 'sc'
        dMdiv = 0.8576*Cl^3 - 1.8139*Cl^2 + 1.0609*Cl - 0.1775;
    end
    
    Mdiv = Mcruise + 0.004 - dMdiv;
    tc = airfoilSweepCurve(Mdiv); %thickness to chord
    c = cosd(sweep)^2*(tc^2)*AR;
    Cl_L = CLestimateLand(c); %low speed clmax during landing
    Cl_TO = CLestimateTO(c); %low speed clmax during takeoff
    WS_L = (Vao/1.3)^2 * sigma*Cl_L/296; %wing loading during approach
    
    Vcruise = Mcruise * 576.4; %cruising speed
    Rao = range + 200 + 0.75*Vcruise; %Range all out
    f_JT8D = 5.1914798324E-13*Rao^3 - 1.1763321777E-08*Rao^2 + 1.2398883131E-04*Rao + 7.6679432694E-03; %fuel fraction for JT8D engines
    f_JT9D = f_JT8D * 0.61/0.78; %fuel fraction for JT9D engines
    f_JT9D = f_JT9D * 1.04; %ratio to adapt
    %f_JT9D = f_JT8D; %if JT8D engines are used
    X = 1; %fuel consumed
    WS_TO = WS_L/(1-X*f_JT9D); %wing loading takeoff
    WS_IC = WS_TO * 0.965; %wing loading initial cruise
    Cl_IC = WS_IC/(1481*PR*Mcruise^2) %cl initial cruise

