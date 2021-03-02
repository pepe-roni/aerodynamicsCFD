%% Spanwise Geometry Plotting
clear str

str = 'Y';

if str == 'Y'
    
    AFoil = importdata('sources/NACA4415_Geometry.dat'); % Import Airfoil (must be in Selig format with no header)
    AFoil(:,1) = -AFoil(:,1); % Flip airfoil X for better plotting
    
    figure(2)
    for i = 1:length(radialStation)
        rotMatrixSpan = [cos(deg2rad(beta(i))) sin(deg2rad(beta(i))); % Rotation matrix blade twist
            -sin(deg2rad(beta(i))) cos(deg2rad(beta(i)))];
        
        localAFoilSection = chord(i)*AFoil;
        localAFoilSection(:,1) = localAFoilSection(:,1) + ones(length(AFoil),1)*.25*chord(i);
        localSection = localAFoilSection*rotMatrixSpan;
        plot(localSection(:,1),localSection(:,2),'Color',[0, .1, i/length(radialStation)]);
        hold on
        
    end
    
    grid on
    grid minor
    hold off
else
end

hold off