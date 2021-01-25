clear 
clc
close all

%freestream velocity
Vinf = 1;
%panels
n =8;

%std geometry 
x = zeros(1,360);
y = zeros(1,360);
for i=1:360
    x(i) = sind(i);
    y(i) = cosd(i);
end
figure(1)
plot(x,y)
hold on

%find verticies
xVertex = [];
yVertex = [];
initAngle = 360/(2*n);
for i=initAngle:(360/n):(360+(360/n))
    xVertex(end+1) = sind(i);
    yVertex(end+1) = cosd(i);
end
plot(xVertex,yVertex)

%solve midpoints 
xMidpoint = zeros(1,n);
yMidpoint = zeros(1,n);
for i=1:numel(xVertex)-1
    xMidpoint(i) = ((xVertex(i)+xVertex(i+1))/2);
    yMidpoint(i) = ((yVertex(i)+yVertex(i+1))/2);
end
plot(xMidpoint,yMidpoint,'o')
legend("Actual Body Geometry","Source Panel Distribution","Panel Midpoints")
title("Geometry of Body and Source Panel n="+n)
xticks(-1:0.2:1)
yticks(-1:0.2:1)
xlim([-1.1 1.1])
ylim([-1.1 1.1])

%find phi and B
phi = zeros(1,n);
for i=1:numel(xMidpoint)
    beta(i)= atan2(yMidpoint(i),xMidpoint(i));
    phi(i) = beta(i)-(pi/2);
end

%calculations for integral matrix (nxn)
integrals = zeros(n);
integrals2 = zeros(n);

for i=1:n
    for j=1:n
        if i==j
            integrals(i,j) = pi;
            integrals2(i,j) = 0;
        else
            A = -(xMidpoint(i)-xVertex(j))*cos(phi(j))-(yMidpoint(i)-yVertex(j))*sin(phi(j));
            B = (xMidpoint(i)-xVertex(j))^2 + (yMidpoint(i)-yVertex(j))^2;
            C = sin(phi(i)-phi(j));
            D = (yMidpoint(i)-yVertex(j))*cos(phi(i))-(xMidpoint(i)-xVertex(j))*sin(phi(i));
            Sj = sqrt((xVertex(j+1)-xVertex(j))^2 + (yVertex(j+1)-yVertex(j))^2);
            E = sqrt(B-A^2);
            integrals(i,j) = (C/2)*log((Sj^2+2*A*Sj+B)/B)+((D-A*C)/E)*(atan2((Sj+A),E)-atan2(A,E));
            integrals2(i,j) = ((D-A*C)/(2*E))*log((Sj^2+2*A*Sj+B)/B)-C*(atan2(Sj+A,E)-atan2(A,E));
        end
    end
end

%solve for free stream velocity components
Vinfs = zeros(n,1);
Vi = zeros(n,1);
Vn = zeros(n,1);
for i=1:n
    Vinfs(i) = Vinf*sin(beta(i));
    Vn(i) = Vinf*cos(beta(i));   
end

%source strengths cancel, solve using Ax=b, x = A\b 
%where A is the coefficient matrix, I(i,j)/2pi and b is -Vn
%3.150(add half lambda to j=i)
coefficient = integrals/(pi*2);
for i=1:n
    coefficient(i,i) = 0.5;
end
%solving Ax=b where x is lambda
lambda = coefficient\(-Vn);

%check 3.157(source strengths = 0); Sj is equal for all panels here
Sj = zeros(n,1); j=1;
Sj(:,1) = sqrt((xVertex(j+1)-xVertex(j))^2 + (yVertex(j+1)-yVertex(j))^2); 
str = dot(lambda,Sj);
disp("Sum of source strengths is "+ str)%very close to 0

%solve for velocity 3.156
Vi = Vinfs + integrals2*lambda/(2*pi); 

%solve for cp 3.38
Cp = 1-(Vi/Vinf).^2;

hold off
figure(2)
hold on
plot(beta+pi,Cp,"o") %we will shift 180 deg 

%predicted cp 3.101
xPred=[0:0.025:2*pi];
yPred=1-4*sin(xPred).^2;
plot(xPred,yPred)

title("Pressure Coefficient of Non-lifting Cylinder: SPM(n="+n+") vs Analytical")
ylabel("Pressure Coefficient Cp")
xlabel("Radians")
legend("Source Panel Mesh","Analytical Calculations")
xlim([0 2*pi])
ylim([-3.5 1.5])
hold off


