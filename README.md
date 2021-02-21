# aerodynamicsCFD
Matlab scripts to calculate, analyze, and simulate aerodynamic flows of various sorts. 

# 195 Propeller Analysis
In this course, I analyzed propeller geometry and its impact on aerodynamic properties. It follows the "Design of Optimum Propellers" from AIAA-83-0190, and I was taught by none other than the author himself, Robert Liebeck.

propellerAnalysis.m analyzes the NACA 4415 airfoil with given properties of onset velocity, diameter, RPM, and num of blades. The characteristics of the airfoil were also digitized and interpolated in the MATLAB script, which is able to get CL values from angle of attack and CD from CL. The script utilizes iterative convergence to calculate flow angle (phi), axial interference factor (A/a1), and rotational interference factor (AP/a2). General aerodynamic properties are also calculated, such as CL, L/D, Reynolds Number, and Mach Number. 

These results were then used to numerically integrate at each blade station. 

<a href="https://www.codecogs.com/eqnedit.php?latex=C_{T}=\frac{1}{\rho&space;n&space;^2&space;D^4}\int_{R_{0}}^{R}&space;\frac{1}{2}&space;\rho&space;W^2&space;BcC_{y}\delta&space;r" target="_blank"><img src="https://latex.codecogs.com/gif.latex?C_{T}=\frac{1}{\rho&space;n&space;^2&space;D^4}\int_{R_{0}}^{R}&space;\frac{1}{2}&space;\rho&space;W^2&space;BcC_{y}\delta&space;r" title="C_{T}=\frac{1}{\rho n ^2 D^4}\int_{R_{0}}^{R} \frac{1}{2} \rho W^2 BcC_{y}\delta r" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=C_{T}=\frac{1}{\rho&space;n&space;^2&space;D^4}\int_{R_{0}}^{R}&space;\frac{1}{2}&space;\rho&space;W^2&space;BcC_{y}\delta&space;r{\color{Orange}&space;}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?C_{T}=\frac{1}{\rho&space;n&space;^2&space;D^4}\int_{R_{0}}^{R}&space;\frac{1}{2}&space;\rho&space;W^2&space;BcC_{y}\delta&space;r{\color{Orange}&space;}" title="C_{T}=\frac{1}{\rho n ^2 D^4}\int_{R_{0}}^{R} \frac{1}{2} \rho W^2 BcC_{y}\delta r{\color{Orange} }" /></a>
