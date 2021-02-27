clear 
clc
close all
warning off

%init conditions
rho = 23.77*10^-4; %density of air at sea level
%rho = 12.67*10^-4; %density 20000ft
%rho = 0.00149620;
mu = 3.737*10^-7; %viscosity
%mu = 0.00022927;

%%DEFINING VARIABLES 
D = 5.75;
R = D/2;
B = 2;
v = 95;
rpm = 2400;
thrust = 207;

% D = 14;
% R = D/2;
% B = 4;
% v = 400*88/60; %ft/s
% rpm = 1200;
% thrust = 3435.4;

n = rpm/60;
omega = rpm*pi/30;


%load external geometry for propeller, change this to analyze other geo
propellerGeo = readtable('propellerGeometry.csv');
%propellerGeo = readtable('propellerGeometryDesignRR.csv');
%propellerGeo = readtable('propellerGeometryP51.csv');
radius = propellerGeo.Root_ft_;
chord = propellerGeo.Chord_ft_;
beta = propellerGeo.Beta_deg_;

%load external lift drag AOT data (if needed) EQs defined @bottom
perfLift = readtable('performanceLift2.csv');
global alpha_data
global c_lift
global c_lift2
global c_drag
alpha_data = perfLift.ALPHA;
c_lift = perfLift.CL;
perfDrag = readtable('performanceDrag2.csv');
c_lift2 = perfDrag.CL;
c_drag = perfDrag.CD;

fprintf(' I    R     CHORD    BETA    PHI      CCL    L/D      RN          MACH   A        AP\n\n');

%%analysis procedure

for j = 1:125
    %reset init values
    dBeta = 0;
    thrust_i = 0;
    phi1 = 0;
    phi2 = 5; 
    a1 = 0; 
    a2 = 0;
    error = 0.0001;
    WW = zeros(21,1);
    CY = zeros(21,1);
    CX = zeros(21,1);
    CL = zeros(21,1);
    CD = zeros(21,1);
    S = zeros(21,1);
    xi = zeros(21,1);
    ALPHA  = zeros(21,1);
%     while thrust_i<thrust
    clc
    fprintf(' I    R     CHORD    BETA    PHI      CCL    L/D      RN          MACH   A        AP\n\n');
    for i=1:numel(beta)
        phi1 = atan2((v*(1+a1)),(omega*radius(i)*(1-a2)));
        while abs(phi2 - phi1)>error
            beta(i) = deg2rad(beta(i)+dBeta);
            alpha = beta(i) - phi1;
            alphaDeg = rad2deg(alpha);
            W = v*(1+a1)/sin(phi1);
            Re = (rho*W*chord(i))/mu;  
            cy = alphaLift(alphaDeg)*cos(phi1)-liftDrag(alphaLift(alphaDeg))*sin(phi1); %thrust force coeff
            cx = alphaLift(alphaDeg)*sin(phi1)+liftDrag(alphaLift(alphaDeg))*cos(phi1); %torque force coeff

            %tip correction
            if radius(i)/R == 1
                radius(i) = radius(i) - 1e-15;
            end

            phiT = atan2(tan(phi1)*radius(i),R);
            f = (B/2)*((1-radius(i)/R)/sin(phiT));
            F = (2/pi)*acos(exp(-f)); %prant loss 

            sigma = B*chord(i)/(2*pi*radius(i));
            a1 = (sigma/(4*F))*(cy/sin(phi1)^2)/(1-(sigma/(F*4))*(cy/sin(phi1)^2));
            a2 = ((sigma/(4*F))*(cx/(sin(phi1)*cos(phi1))))/((1+(sigma/(F*4))*(cx/(sin(phi1)*cos(phi1)))));

            if abs(a1) > .7 || abs(a2) > .7
                a2 = .4;
            elseif isnan(a1) || isnan(a2)
                f = .001;
                F = .001;
            end
            phi2 = phi1;
            phi1 = atan2((v*(1+a1)),(omega*radius(i)*(1-a2)));
            phi1 = phi2 + 0.4*(phi1-phi2);
            phiDeg = rad2deg(phi1);
            beta(i) = rad2deg(beta(i));
            mach = W/1100;
        end


        fprintf('%2i   %.2f   %.4f   %.2f   %.3f   %.2f   %.2f',i, radius(i), chord(i), beta(i), phiDeg, alphaLift(alphaDeg), alphaLift(alphaDeg)/liftDrag(alphaLift(alphaDeg)));
        fprintf('   %10.2f   %.2f   %.4f   %.4f\n', Re, mach, a1, a2)

        WW(i,1) = W;
        CY(i,1) = cy;
        CX(i,1) = cx;
        CL(i,1) = alphaLift(alphaDeg);
        CD(i,1) = liftDrag(alphaLift(alphaDeg));
        S(i,1) = B*chord(i)/(pi*R^2);
        xi(i,1) = radius(i)/R;
        ALPHA(i,1) = alphaDeg;

    end

    dR = radius;
    integrand_T = 0.5*rho*WW.^2*B.*chord.*CY;
    thrust_i = trapz(dR, integrand_T); %lbs?
    ct = thrust_i/(rho*n^2*D^4); 
    integrand_P = 0.5*omega*rho*WW.^2*B.*chord.*CX.*radius;
    power = trapz(dR, integrand_P); %this is in ft lbs/sec
    cp = power/(rho*n^3*D^5);
    hp = power/550; %convert from ft lbs/sec to hp
    AR = v/(n*D); %advance ratio, J
    ETA = ct*AR/cp;
    solidity = trapz(dR,S);
    AF = (100000/(16*D))*trapz(xi,chord.*xi.^3); %activity factor for a single blade
    AF = AF*B; %activity factor for the entire propeller
   
%     dBeta = dBeta+ 0.0001;
%     end
    CP(j) = cp;
    CT(j) = ct;
    J(j) = v/(n*D); %advance ratio, J
    eta(j) = ct*J(j)/cp;
    v = v+1;
end
fprintf('\nThrust: %.2f  CT: %.4f  Power: %.1f  CP: %.4f  HP: %.2f  AdvR: %.3f  ETA: %.4f\n',thrust_i,ct,power,cp,hp,AR,ETA);
fprintf('Solidity: %.3f  AF: %.2f   dBeta: %.5f\n', solidity, AF, dBeta);
%run external functions
plotGeometry
hold off
figure(2)
subplot(2,1,1)
plot(J,eta)

xlabel('J')
ylabel('ETA')
xlim([0.3 1])
ylim([0 1])
title('NACA4415 RPM=2400 V= 95 to 220 ft/s')
subplot(2,1,2)
hold on
plot(J,CP)
plot(J,CT)
xlim([0.3 1])
ylim([0 0.1])
ylabel('CT,CP')
xlabel('J')
legend('Cp','Ct')

warning on

%functions of alpha-lift curve and cl-cd calculated from perf chart
function lift = alphaLift(angle)   
    %lift = -0.0003*angle^3 + 0.0023*angle^2 + 0.1005*angle + 0.3852; %digitize
    %lift = -0.0001*angle^3 - 0.0008*angle^2 + 0.1094*angle + 0.4612;
    %lift = 0.7; %debug
    
    global alpha_data
    global c_lift
    
    %interpolation
    difference = zeros(numel(alpha_data),1);
    for i=1:numel(alpha_data)
        difference(i,1) = abs(alpha_data(i) - angle);
    end
    [~,closestInd] = min(difference);
    difference(closestInd) = inf; %set the value we just found to inf so that we can find 2nd closest
    [~,closestInd2] = min(difference);
    
    lift = interpo(alpha_data(closestInd),alpha_data(closestInd2), angle, c_lift(closestInd), c_lift(closestInd2));
    
end

function drag = liftDrag(l)
    %drag = 0.0097*l^3 - 0.0053*l^2 - 0.004*l + 0.0111; %digitized data
    global c_lift2
    global c_drag
    
    %interpolation
    difference = zeros(numel(c_lift2),1);
    for i=1:numel(c_lift2)
        difference(i,1) = abs(c_lift2(i) - l);
    end
    [~,closestInd] = min(difference);
    difference(closestInd) = inf; %set the value we just found to inf so that we can find 2nd closest
    [~,closestInd2] = min(difference);
    
    drag = interpo(c_lift2(closestInd),c_lift2(closestInd2), l, c_drag(closestInd), c_drag(closestInd2));
end

function interpValue = interpo(x1,x2,x3,y1,y2)
    interpValue = (((x2-x3)*y1) + ((x3-x1)*y2))/(x2-x1);
end