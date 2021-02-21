# aerodynamicsCFD
Matlab scripts to calculate, analyze, and simulate aerodynamic flows of various sorts. 

# 195 Propeller Analysis
In this course, I analyzed propeller geometry and its impact on aerodynamic properties. It follows the "Design of Optimum Propellers" from AIAA-83-0190, and I was taught by none other than one of the authors himself, Robert Liebeck.

propellerAnalysis.m analyzes the NACA 4415 airfoil with given properties of onset velocity, diameter, RPM, and num of blades. The characteristics of the airfoil were also digitized and interpolated in the MATLAB script, which is able to get CL values from angle of attack and CD from CL. The script utilizes iterative convergence to calculate flow angle (phi), axial interference factor (A or a1), and rotational interference factor (AP or a2). General aerodynamic properties are also calculated, such as CL, L/D, Reynolds Number (RN or Re), and Mach Number. 

These results were then used to numerically integrate at each blade station as follows (in orange for both dark and light mode users <3):

<img src="https://latex.codecogs.com/svg.latex?{\color{Orange}&space;C_{T}=\frac{1}{\rho&space;n&space;^2&space;D^4}\int_{R_{0}}^{R}&space;\frac{1}{2}&space;\rho&space;W^2&space;BcC_{y}\delta&space;r}" title="{\color{Orange} C_{T}=\frac{1}{\rho n ^2 D^4}\int_{R_{0}}^{R} \frac{1}{2} \rho W^2 BcC_{y}\delta r}" />

<img src="https://latex.codecogs.com/svg.latex?{\color{Orange}&space;C_{p}=\frac{1}{\rho&space;n&space;^3&space;D^5}\int_{R_{0}}^{R}&space;\frac{1}{2}\Omega\rho&space;W^2&space;BcC_{x}r\delta&space;r&space;}" title="{\color{Orange} C_{p}=\frac{1}{\rho n ^3 D^5}\int_{R_{0}}^{R} \frac{1}{2}\Omega\rho W^2 BcC_{x}r\delta r }" />

Where rho is the density of air, n is propeller rotations per second, omega is propeller radians per second, D is the diameter of the propeller, R is the radius of the propeller tip, R0 is the radius of the propeller at the root, W is the local velocity, B is the number of blades, c is the chord at the specific blade station, Cy is the thrust force coefficent, and Cx is the torque force coefficient.

So the script will then calculate the Thrust (T), Power (P), Coefficent of Thrust (Ct), Coefficient of Power (Cp), Advance Ratio (AR or J), Activity Factor (AF), Efficiency Factor (ETA), and Solidity.

An output can be seen here:
https://github.com/pepe-roni/aerodynamicsCFD/blob/main/195/output.png

# 136 Source Panel Method
The source panel method is useful as a numerical method for solving flows of arbitrary shapes. For this to be possible, the body must be broken into panels of length dS and given a source strength to simulate a solid body; the more panels, the better the resolution of the final result [1]. In this project, the source panel method will be used to analyze flow over a non lifting circular cylinder of radius 1 and compare our results to the proven analytical method. Example 3.17 from Fundamentals of Aerodynamics by John Anderson will be used as a guide.

Increasing the number of panels provided more data points along the analytical calculation plot. It improved the resolution of the plotted points but not the accuracy of the plot. The source panel method proved to match the analytical derivation essentially perfectly. If this was conducted in the real world, the results would likely be similar but not exact because it is difficult to keep the ideal conditions necessary for this method.

A write-up can be seen here: http://sourcepanelmatlab.github.io/
