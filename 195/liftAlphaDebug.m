perfLift = readtable('performanceLift2.csv');
global alpha_data
global c_lift
global c_lift2
global c_drag
alpha_data = perfLift.ALPHA;
c_lift = perfLift.CL;

liftAlpha(0.7)

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

function interpValue = interpo(x1,x2,x3,y1,y2)
    interpValue = (((x2-x3)*y1) + ((x3-x1)*y2))/(x2-x1);
end