# aerodynamicsCFD
Matlab scripts to calculate, analyze, and simulate aerodynamic flows of various sorts. 

# 159 Aircraft Design
<img src="https://github.com/pepe-roni/aerodynamicsCFD/blob/main/MAE%20159/Screenshot%202022-06-15%20143608.jpg" />

The objective of this project was to design a subsonic commercial transport. Composite, hybrid, and aluminum materials were compared as well as design variables such as engine type, engine configuration, aspect ratio, sweep, aisles, seats abreast, cargo weight, and airfoil type. Analysis and preliminary design code was written in MATLAB to form the basis of the aircraft. With texts from Shevell and Schaufele, the mechanical design and drawings of the aircraft was designed in Solidworks and rendered in NX as seen in the report.

With the specifications, three aircraft designs will be presented and analyzed using the design code. Design variables will be changed to observe their impact on performance and risk, ultimately leading to which variables are most indicative of success. For this comparison the baseline will be an all-aluminum airframe and composite-aluminum hybrids and all composite designs will be analyzed. The primary comparison will be factors that affect DOC (direct operating costs) and weight. 

The final design of the aircraft is capable of delivering 4000 lbs of cargo and 200 passengers over 4300 nmi. This can be extended to just below 6000 nmi with the removal of cargo and passenger load. This range exceeds that of the Boeing 737. Composite airframe structures and advancements in turbofan engines make significant improvements to aircraft design. With these additions, the airframe can achieve greater efficiency, lower DOC, and lower fuel and takeoff weights. When coupling this with other options such as high aspect ratio swept wings, supercritical airfoils, winglets, and anti-shock bodies, the improvements add up  significantly. The fully dimensioned plane resembles that of a medium class transport like the 737, and its performance, cost, and other parameters are also extremely similar.

See a full report here! https://github.com/pepe-roni/aerodynamicsCFD/blob/main/MAE%20159/combinepdf%20(1).pdf

# 195 Propeller Analysis
In this project, I analyzed propeller geometry and its impact on aerodynamic properties. It follows the "Design of Optimum Propellers" from AIAA-83-0190, and I was taught by none other than one of the authors himself, Robert Liebeck.

propellerAnalysis.m analyzes the NACA 4415 airfoil with given properties of onset velocity, diameter, RPM, and num of blades. The characteristics of the airfoil were also digitized and interpolated in the MATLAB script, which is able to get CL values from angle of attack and CD from CL. The script utilizes iterative convergence to calculate flow angle (phi), axial interference factor (A or a1), and rotational interference factor (AP or a2). General aerodynamic properties are also calculated, such as CL, L/D, Reynolds Number (RN or Re), and Mach Number. 

These results were then used to numerically integrate at each blade station as follows (in orange for both dark and light mode users <3):

<img src="https://latex.codecogs.com/svg.latex?{\color{Orange}&space;C_{T}=\frac{1}{\rho&space;n&space;^2&space;D^4}\int_{R_{0}}^{R}&space;\frac{1}{2}&space;\rho&space;W^2&space;BcC_{y}\delta&space;r}" title="{\color{Orange} C_{T}=\frac{1}{\rho n ^2 D^4}\int_{R_{0}}^{R} \frac{1}{2} \rho W^2 BcC_{y}\delta r}" />

<img src="https://latex.codecogs.com/svg.latex?{\color{Orange}&space;C_{p}=\frac{1}{\rho&space;n&space;^3&space;D^5}\int_{R_{0}}^{R}&space;\frac{1}{2}\Omega\rho&space;W^2&space;BcC_{x}r\delta&space;r&space;}" title="{\color{Orange} C_{p}=\frac{1}{\rho n ^3 D^5}\int_{R_{0}}^{R} \frac{1}{2}\Omega\rho W^2 BcC_{x}r\delta r }" />

Where rho is the density of air, n is propeller rotations per second, omega is propeller radians per second, D is the diameter of the propeller, R is the radius of the propeller tip, R0 is the radius of the propeller at the root, W is the local velocity, B is the number of blades, c is the chord at the specific blade station, Cy is the thrust force coefficent, and Cx is the torque force coefficient.

So the script will then calculate the Thrust (T), Power (P), Coefficent of Thrust (Ct), Coefficient of Power (Cp), Advance Ratio (AR or J), Activity Factor (AF), Efficiency Factor (ETA), and Solidity.

An output can be seen here:
https://github.com/pepe-roni/aerodynamicsCFD/blob/main/195/output.png

See a full report here! https://docs.google.com/document/d/1Q3HVWSurWS-HKOuec8uOBUdu3PdFcdoMMpoL0xtbG_E/edit?usp=sharing

# 112 Propulsion Systems - Supersonic Turbofan Engine Design
The objective of this project was to design a turbofan engine for a commercial application. In this real life example, we are tuning our parameters to be on par with an example from Boom Technologies (https://boomsupersonic.com). For this project, we will be designing to meet the requirements for their supersonic business jet, the Overture. To have a compelling option for them, we must meet their design specifications which we will cover in the introduction. Additionally, we must tune parameters of the turbofan engine on design characteristics to achieve the desired outcome. These include turbine inlet temperature, compressor pressure ratio, bypass ratio, bypass pressure ratio, and engine inlet diameter. Tradeoffs and compromises must be made in order for all specifications and requirements to be met.

See full report: https://docs.google.com/document/d/1nposa372G0WSS5NIN29Xs4I7jhmCFPcny86xvYq9tRo/edit?usp=sharing

# 112 Propulsion Systems - Ramjet On-design Performance
Ramjets are simple air-breathing engines that operate primarily (and most optimally) in supersonic flight. While subsonic speeds are possible, we will see how the ramjet performs with varying and increasing Mach number flows. The ramjet is largely dependent on the maximum temperature of the combustion chamber. This is typically limited by the constraints and choice of the material. Having a higher temperature in combustion will typically result in better efficiency in both thrust and fuel consumption. In this report, we will explore the thrust, fuel consumption, efficiencies, and other on-design parameters of the ramjet engine in super and partially hypersonic flight.

See full report: https://docs.google.com/document/d/1Ez_rfg51uyH499CT98kBn2j6Q9d20e1HbWmX-NeTIis/edit?usp=sharing

# 136 Compressible Flow - Subsonic and Supersonic Windtunnel Designs
The objective of this project was to design a subsonic wind tunnel with some predefined design requirements. In our case, we have a predetermined test chamber area size, minimum test chamber length, and a required test chamber mach number. The method included analyzing the isentropic flow after the fan section, choosing the areas of the inlet and outlet, and utilizing conservation of momentum and energy. Factors to keep in mind while designing are size, cost, and manufacturability considerations, aside from the actual engineering factors required. For instance, a powerful fan may be inefficient due to losses and consume excessive energy while having a very high cost to manufacture. The results of our design was a reasonably sized wind tunnel capable of the mentioned design requirements while performing relatively efficiently due to optimal cross sectional areas and therefore fan choice. 

See Subsonic full report: https://docs.google.com/document/d/16VsFLv-CrFGwO4kMxLl0i2oblSZTRtRvYnd__ujW6EA/edit?usp=sharing

See Supersonic full report: https://docs.google.com/document/d/15yiEJpc4Wh3ik4uxX1nBZMRepogn_o8ck7cygxznLxI/edit?usp=sharing


# 136 Source Panel Method
The source panel method is useful as a numerical method for solving flows of arbitrary shapes. For this to be possible, the body must be broken into panels of length dS and given a source strength to simulate a solid body; the more panels, the better the resolution of the final result [1]. In this project, the source panel method will be used to analyze flow over a non lifting circular cylinder of radius 1 and compare our results to the proven analytical method. Example 3.17 from Fundamentals of Aerodynamics by John Anderson will be used as a guide.

Increasing the number of panels provided more data points along the analytical calculation plot. It improved the resolution of the plotted points but not the accuracy of the plot. The source panel method proved to match the analytical derivation essentially perfectly. If this was conducted in the real world, the results would likely be similar but not exact because it is difficult to keep the ideal conditions necessary for this method.

A write-up can be seen here: https://docs.google.com/document/d/1TA7tnR35MIVRT1I85CGtVVVetPy_QTQMVCb3TI2a4Mk/edit?usp=sharing
