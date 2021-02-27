clear 
clc
close all
warning off

%init conditions
rho = 23.77*10^-4; %density of air at sea level
rho = 12.67*10^-4; %density 20000ft
mu = 3.737*10^-7; %viscosity
%mu = 0.00022927;

%%DEFINING VARIABLES 
D = 14;
R = D/2;
B = 4;
v = 400*88/60; %ft/s
pwr = 4000; %bhp
pwr = pwr*550; %lbft/s
thrust = 0;
rpm = 1200;
rroot = 1.5/2; %radius of the root

Cp = 2*pwr/(rho*v^3*pi*R^2);
Ct = 2*thrust/(rho*v^2*pi*R^2);

n = rpm/60;
omega = rpm*pi/30;

lambda = v/(omega*R); %speed ratio

%load external geometry for propeller, change this to analyze other geo
propellerGeo = readtable('propellerGeometry.csv');
%propellerGeo = readtable('propellerGeometryP51.csv');
radius(:,1) = (rroot:(R-rroot)/20:R);
% chord = propellerGeo.Chord_ft_;
% beta = propellerGeo.Beta_deg_;
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

%%Design
zeta1 = 0;
zeta2 = 5; 
a1 = 0; 
a2 = 0;
error = 0.0001;

xi = zeros(21,1);
X = zeros(21,1);
phi = zeros(21,1);
F = zeros(21,1);
alphaRad = zeros(21,1);
beta = zeros(21,1);
chord = zeros(21,1);
I1P = zeros(21,1);
I2P = zeros(21,1);
J1P = zeros(21,1);
J2P = zeros(21,1);
while abs(zeta2 - zeta1)>error
    
    for i=1:numel(radius)
         %tip correction
        if radius(i)/R == 1
            radius(i) = radius(i) - 1e-15;
        end
        
        xi(i) = radius(i)/R;
        X(i) = omega*radius(i)/v;
        phiT = atan2(lambda*(1+zeta1/2),1);
        phi(i)= atan2(tan(phiT),xi(i)); 
        f = (B/2)*((1-radius(i)/R)/sin(phiT));
        F(i,1) = (2/pi)*acos(exp(-f));
        G(i,1) = F(i)*X(i)*sin(phi(i))*cos(phi(i));
        
        Cl = 0.7; %design criteria
        alpha = liftAlpha(0.7);
        Cd = liftDrag(Cl);
        e = Cd/Cl;
        alphaRad(i) = deg2rad(alpha);
        
        gamma = 2*pi*v^2*zeta1*G/(R*omega);
        wt=v*zeta1*sin(phi(i))*cos(phi(i));
        wn=wt/sin(phi(i));
        Wc(i,1) = (4*pi*v^2*G(i)*zeta1)/(B*omega*Cl);
        
        %interference factors
        a1 = ((zeta1/2)*cos(phi(i))^2)*(1-e*tan(phi(i)));
        a2 = (zeta1/(2*X(i)))*cos(phi(i))*sin(phi(i))*(1-e/tan(phi(i)));
        W = v*(1+a1)/sin(phi(i));
        
        %geometry
        beta(i,1) = alphaRad(i)+phi(i);
        chord(i,1) = Wc(i)/W;
        
        I1P(i,1) = 4*xi(i)*G(i)*(1-e*tan(phi(i)));
        I2P(i,1) = lambda*(I1P(i)/2*xi(i))*(1+e/tan(phi(i)))*sin(phi(i))*cos(phi(i));
        J1P(i,1) = 4*xi(i)*G(i)*(1+e/tan(phi(i)));
        J2P(i,1) = (J1P(i)/2)*(1-e*tan(phi(i)))*cos(phi(i))^2;
    end
   zeta2 = zeta1;
   I1 = trapz(xi,I1P);
   I2 = trapz(xi,I2P);
   J1 = trapz(xi,J1P);
   J2 = trapz(xi,J2P);
   
   if Cp == 0 %this is if thrust was specified
       zeta1=I1/(2*I2)-(((I1/(2*I2))^2)-(Ct/I2))^(1/2);
       Cp = J1*zeta1 + J2*(zeta1^2);
   else %this is if power was specified
       zeta1 = -(J1/(2*J2))+(((J1/(2*J2))^2 +(Cp/J2)))^(1/2);
       Ct = I1*zeta1 - I2*(zeta1^2);
   end
end
beta = rad2deg(beta);
fprintf('\nDesign Conditions:\n D:%ift  B:%i  V:%4.1fft/s  RPM:%i   rho:%.2ilbf/ft^3  Desired CL:%.3f\n\n',D,B,v,rpm,rho,Cl)
fprintf(' I   Radius   Chord     Beta\n')
for i=1:numel(chord)
    fprintf('%2i   %.4f   %.4f    %.4f\n',i,radius(i),chord(i),beta(i))
end
thrust = Ct*(rho*v^2*R^2*pi*0.5);
power = Cp*0.5*rho*v^3*pi*R^2;
power = power/550;
eta = Ct/Cp;
AR = v/(n*D); %advance ratio, J
AF = (100000/(16*D))*trapz(xi,chord.*xi.^3); %activity factor for a single blade
AF = AF*B;

fprintf('\nCt:%.3f  Thrust:%.1f  Cp:%.3f   Power(HP):%.1f   ETA:%.2f   AR:%.2f    AF:%.2f\n',Ct,thrust,Cp,power,eta,AR,AF)

%run external functions
plotGeometry

warning on

%functions of alpha-lift curve and cl-cd calculated from perf chart
function alpha = liftAlpha(lift)   
    %lift = -0.0003*angle^3 + 0.0023*angle^2 + 0.1005*angle + 0.3852; %digitize
    %lift = -0.0001*angle^3 - 0.0008*angle^2 + 0.1094*angle + 0.4612;
    %lift = 0.7; %debug
    
    global alpha_data
    global c_lift
    
    %interpolation
    difference = zeros(numel(c_lift),1);
    for i=1:numel(c_lift)
        difference(i,1) = abs(c_lift(i) - lift);
    end
    [~,closestInd] = min(difference);
    difference(closestInd) = inf; %set the value we just found to inf so that we can find 2nd closest
    [~,closestInd2] = min(difference);
    
    alpha = interpo(c_lift(closestInd),c_lift(closestInd2), lift, alpha_data(closestInd), alpha_data(closestInd2));
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