% Paste me at the bottom of the code. I will plot axial geometry, spanwise (blade stations), and 3D geometry.
%The 3D geometry is a bit experimental so you might have to adjust some thing accordingly!

%% Geometry Plotting
clear str
close all

str = 'Y';

blades = B;
radialStation = radius;
hubRadius = radius(1);
if str == 'Y'
    bladeSpace = 2*pi/blades; % Blade spacing
    
    if mod(blades,2)~=0 % This checks for odd number of blades
        plotRot = pi/2; % Additional rotation for odd blades
    else
        plotRot = 0;
    end
    figure(1)
    for i = 1:blades
        
        rotMatrix = [cos(bladeSpace*(i-1)- plotRot) sin(bladeSpace*(i-1)- plotRot); % Rotation matrix
            -sin(bladeSpace*(i-1)- plotRot) cos(bladeSpace*(i-1)- plotRot)];
        
        for j = 1:length(radialStation)
            bladePlotLE(1,j) = radialStation(j);
            bladePlotTE(1,j) = radialStation(j);
            bladePlotLE(2,j) = .25*chord(j); % Blade Leading Edge coordinates
            bladePlotTE(2,j) = bladePlotLE(2,j) - chord(j); % Blade Trailing Edge coordinates
            
            bladePlotLE(:,j) = rotMatrix*bladePlotLE(:,j); % Rotate blade plotting
            bladePlotTE(:,j) = rotMatrix*bladePlotTE(:,j); % Rotate blade plotting
            
        end
        plot(bladePlotLE(1,:),bladePlotLE(2,:),'Color','k'); % Plot LE
        hold on
        plot(bladePlotTE(1,:),bladePlotTE(2,:),'Color','r'); % Plot TE
        hold on
        
        
    end
    xlim([-R R])
    ylim([-R R])
    
    viscircles([0,0],hubRadius,'Color','k','LineWidth',.75); % Plots hub
    hold on
    viscircles([0,0],R,'Color','b','LineStyle','--','LineWidth',.5); % Plot tip arc
    grid on
    grid minor
    
else
end