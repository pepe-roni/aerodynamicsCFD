clear str

str = 'Y';
blades = B;

if str == 'Y'
    
    bladeSpace = 2*pi/blades; % Blade spacing
    
    AFoil = importdata('sources/NACA4415_Geometry.dat'); % Import Airfoil (must be in Selig format with no header)
    AFoil(:,1) = -AFoil(:,1); % Flip airfoil X for better plotting
    
    if mod(blades,2)~=0 % This checks for odd number of blades
        plotRot = pi/2; % Additional rotation for odd blades
    else
        plotRot = 0;
    end
    
    figure(3)
    for i = 1:blades
        
        rotMatrixAxial = [cos(bladeSpace*(i-1)- plotRot) 0 sin(bladeSpace*(i-1)- plotRot); % Rotation matrix axial
            0 1 0;
            -sin(bladeSpace*(i-1)- plotRot) 0 cos(bladeSpace*(i-1)- plotRot)];
        
        for i = 1:length(radialStation)
            rotMatrixSpan = [cos(deg2rad(beta(i))) sin(deg2rad(beta(i))); % Rotation matrix blade twist
                -sin(deg2rad(beta(i))) cos(deg2rad(beta(i)))];
            
            localAFoilSection = chord(i)*AFoil;
            localAFoilSection(:,1) = localAFoilSection(:,1) + ones(length(AFoil),1)*.25*chord(i);
            localSection = horzcat(localAFoilSection*rotMatrixSpan,-radialStation(i)*ones(length(AFoil),1))*rotMatrixAxial;
            
            plot3(localSection(:,1),localSection(:,3),localSection(:,2),'Color','k');
            hold on
            
        end
        
        
    end
    
    xlim([-R R])
    ylim([-R R])
    zlim([-R R])
    grid on
    grid minor
    
    % Plot Spinner Cone
    xx = 0:.01:hubRadius;
    spinnerHeight = .4*radius; % "Tuned" to make the spinner look proportional
    yy = sqrt(spinnerHeight*(1-(xx./hubRadius).^2));
    spinnerCone = transpose([xx;yy - yy(end) - chord(1);zeros(1,length(xx))]);
    numConeStations = 30;
    for k = 1:numConeStations
        angle = 2*pi*k/numConeStations;
        rotMatrixSpinner = [cos(angle) 0 sin(angle); % Rotation matrix spinner cone
            0 1 0;
            -sin(angle) 0 cos(angle)];
        drawSpinner = spinnerCone*rotMatrixSpinner;
        plot3(drawSpinner(:,1),drawSpinner(:,3),drawSpinner(:,2),'Color','b')
        hold on
    end
    
else
end