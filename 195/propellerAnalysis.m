clear 
clc
close all

%init conditions
rho = 23.77*10^-4; %density of air at sea level
mu = 1.57*10^-4; %viscosity

%%variables given on page 40
v = 161.33;
rpm = 2400;
B = 2;
omega = rpm*pi/30;
R=5.75/2;

%load external geometry for propeller
propellerGeo = readtable('propellerGeometry.csv');
radius = propellerGeo.Root_ft_;
chord = propellerGeo.Chord_ft_;
beta = propellerGeo.Beta_deg_;
%load external lift drag AOT data (if needed) EQs defined @bottom
perfLift = readtable('performanceLift.csv');
alpha = perfLift.alpha;
c_lift = perfLift.Cl;
perfDrag = readtable('performanceDrag.csv');
c_lift2 = perfDrag.cl;
c_drag = perfDrag.cd;

%%analysis procedure

%functions of alpha-lift curve and cl-cd calculated from perf chart
function lift = alphaLift(angle)
    lift = -0.0003*angle^3 + 0.0023*angle^2 + 0.1005*angle + 0.3852;
end

function drag = liftDrag(l)
    drag = 0.0097*l^3 - 0.0053*l^2 - 0.004*l + 0.0111;
end