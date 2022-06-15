clc
clear
close all

%2022 specs
range = 4000; %range in miles
Vao = 135; %landing approach speed in knots
Mcruise = 0.82; %cruise mach number
airfoilType = 'sc'; 
sweep = 25;
AR = 8;
sigma = 0.953; %density ratio
PR = 0.2360; %pressure ratio
wingEngines = 2;
fuselageEngines = 0;
engineType = 'JT9D';
engineAdv = true;
composite = false;
hybrid = false;
TOFL = 6000; %take off field length in ft
X = 1; %fuel consumed
PAX = 200;
W_cargo = 4000;
abreast = 5;
aisle = 1;
alt = 35000;

%prof specs
% PAX = 275;
% W_cargo = 12000;
% range = 6000;
% TOFL = 9000;
% Vao = 140;
% X = 0.75;

%modified specs
% PAX = 135;
% aisle = 1;
% abreast = 6;
% AR = 7;
% sweep = 32;
% TOFL = 8000;
% range = 2550;
% W_cargo = 3000;
% Mcruise = 0.8;
% Vao = 140;
% fuselageEngines = 3;
% wingEngines = 0;
% engineType = 'JT9D';
% X = 0.75;



for i=1:6
    for j=1:6
        engineCount = wingEngines+fuselageEngines; 
        %1 or 0, use interpolation when using non 25, 30, 35, 40 angles
        if sweep ~= 25 && sweep ~=30 && sweep ~=35 && sweep ~=40
            interpolate = 1;
        else
            interpolate = 0;
        end

        if engineCount == 2
            syms1 = @(x) -2.75700E-01*x^2 + 3.23565E+01*x - 2.17775E+01;
        elseif engineCount == 3
            syms1 = @(x) -3.14122E-01*x^2 + 3.59436E+01*x - 2.10352E+01;
        elseif engineCount == 4
            syms1 = @(x) -3.50586E-01*x^2 + 3.73883E+01*x - 1.29689E+01;
        end

        if interpolate == 1
            if airfoilType == 'c'
                if sweep>25 && sweep<30
                    airfoilSweepCurve = @(x) 1.80081E-02*x^2 - 5.59031E-01*x + 5.28813E-01; %conventional 25 deg
                    airfoilSweepCurve2 = @(x) 1.44744E-02*x^2 - 5.25458E-01*x + 5.12939E-01; %conventional 30 deg
                    sweepValues = [25,30];
                elseif sweep>30 && sweep<35
                    airfoilSweepCurve = @(x) 1.44744E-02*x^2 - 5.25458E-01*x + 5.12939E-01; %conventional 30 deg
                    airfoilSweepCurve2 = @(x) 7.10020E-03*x^2 - 4.78583E-01*x + 4.89769E-01; %conventional 35 deg
                    sweepValues = [30,35];
                elseif sweep>35 && sweep<40
                    airfoilSweepCurve = @(x) 7.10020E-03*x^2 - 4.78583E-01*x + 4.89769E-01; %conventional 35 deg
                    airfoilSweepCurve2 = @(x) -1.58922E-03*x^2 - 4.26203E-01*x + 4.63453E-01; %conventional 40 deg
                    sweepValues = [35,40];
                end
            elseif airfoilType == 'sc'
                if sweep>25 && sweep<30
                    airfoilSweepCurve = @(x) 3.73842E+00*x^2 - 6.73394E+00*x + 3.10823E+00 %supercritical 25 deg
                    airfoilSweepCurve2 = @(x) -1.75962E+01*x^3 + 4.71604E+01*x^2 - 4.26422E+01*x + 1.30762E+01; %supercritical 30 deg
                    sweepValues = [25,30];
                elseif sweep>30 && sweep<35
                    airfoilSweepCurve = @(x) -1.75962E+01*x^3 + 4.71604E+01*x^2 - 4.26422E+01*x + 1.30762E+01; %supercritical 30 deg
                    airfoilSweepCurve2 = @(x) -3.76659E+01*x^3 + 1.00937E+02*x^2 - 9.08776E+01*x + 2.75695E+01; %supercritical 35 deg
                    sweepValues = [30,35];
                elseif sweep>35 && sweep<40
                    airfoilSweepCurve = @(x) -3.76659E+01*x^3 + 1.00937E+02*x^2 - 9.08776E+01*x + 2.75695E+01; %supercritical 35 deg
                    airfoilSweepCurve2 = @(x) -4.46244E+02*x^3 + 1.17289E+03*x^2 - 1.02924E+03*x + 3.01674E+02; %supercritical 40 deg
                    sweepValues = [35,40];
                end
            end
        elseif interpolate == 0
            if airfoilType == 'c'
                if sweep == 25
                    airfoilSweepCurve = @(x) 1.80081E-02*x^2 - 5.59031E-01*x + 5.28813E-01; %conventional 25 deg
                elseif sweep == 30
                    airfoilSweepCurve = @(x) 1.44744E-02*x^2 - 5.25458E-01*x + 5.12939E-01; %conventional 30 deg
                elseif sweep == 35
                    airfoilSweepCurve = @(x) 7.10020E-03*x^2 - 4.78583E-01*x + 4.89769E-01; %conventional 35 deg
                elseif sweep == 40
                    airfoilSweepCurve = @(x) -1.58922E-03*x^2 - 4.26203E-01*x + 4.63453E-01; %conventional 40 deg
                end
            elseif airfoilType == 'sc'
                if sweep == 25
                    airfoilSweepCurve = @(x) 3.73842E+00*x^2 - 6.73394E+00*x + 3.10823E+00; %supercritical 25 deg
                elseif sweep == 30
                    airfoilSweepCurve = @(x) -1.75962E+01*x^3 + 4.71604E+01*x^2 - 4.26422E+01*x + 1.30762E+01; %supercritical 30 deg
                elseif sweep == 35
                    airfoilSweepCurve = @(x) -3.76659E+01*x^3 + 1.00937E+02*x^2 - 9.08776E+01*x + 2.75695E+01; %supercritical 35 deg
                elseif sweep == 40
                    airfoilSweepCurve = @(x) -4.46244E+02*x^3 + 1.17289E+03*x^2 - 1.02924E+03*x + 3.01674E+02; %supercritical 40 deg
                end
            end
        end

        %clmax estimates
        CLestimateLand = @(x) 109.223789*x^3 - 67.377088*x^2 + 16.449711*x + 2.002030; %landing Cl from "c"
        CLestimateTO = @(x) 84.325414*x^3 - 65.343054*x^2 + 17.178855*x + 1.029074;  %TO Cl from "c"

        %shevell drag, aircraft performance
        RNestimateCf = @(x) (7.7538691E+01)*x^(-1.9414863E-01); %approx Cf from RN, only from 10^6 to 10^8
        bodyForm = @(x)  -1.51985E-03*x^3 + 3.89461E-02*x^2 - 3.59054E-01*x + 2.30547E+00; %approx body form from fineness ratio

        %tofl and K relations


        addF = 0; %added fuel
        errRange = 1000;
        err = 1; Cl = 0.5; %inital guess
        addW = 0; %added weight  = 0
        errorWeight = -10;

        while(errorWeight<0)
            while(errRange>5)
                %CL CONVERGENCE
                while(err>0.00011)
                    if airfoilType == 'c'
                        dMdiv = -0.2044*Cl^2 - 0.1122*Cl + 0.1236;
                    elseif airfoilType == 'sc'
                        dMdiv = 0.8576*Cl^3 - 1.8139*Cl^2 + 1.0609*Cl - 0.1775;
                    end

                    Mdiv = Mcruise + 0.004 - dMdiv;

                    if interpolate == 1
                        tc1 = airfoilSweepCurve(Mdiv);
                        tc2 = airfoilSweepCurve2(Mdiv);
                        tc = tc1*((sweepValues(2)-sweep)/(sweepValues(2)-sweepValues(1))) + tc2*(1 - ((sweepValues(2)-sweep)/(sweepValues(2)-sweepValues(1))));
                    elseif interpolate == 0
                        tc = airfoilSweepCurve(Mdiv); %thickness to chord
                    end

                    c = cosd(sweep)^2*(tc^2)*AR;
                    Cl_L = CLestimateLand(c); %low speed clmax during landing
                    Cl_TO = CLestimateTO(c); %low speed clmax during takeoff
                    WS_L = (Vao/1.3)^2 * sigma*Cl_L/296; %wing loading during approach

                    Vcruise = Mcruise * 576.4; %cruising speed
                    Rao = range + 200 + 0.75*Vcruise; %Range all out
                    f_JT8D = 5.1914798324E-13*Rao^3 - 1.1763321777E-08*Rao^2 + 1.2398883131E-04*Rao + 7.6679432694E-03; %fuel fraction for JT8D engines

                    if engineType == 'JT8D'
                        f_JT9D = f_JT8D + addF; %if JT8D engines are used
                        T_SLS = 14500;
                        T_m = @(x) (14.47*x^4 - 41.021*x^3 + 33.863*x^2 - 14.296*x + 14.646)*1000;
                    elseif engineType == 'JT9D'
                        f_JT9D = f_JT8D * 0.61/0.78; %fuel fraction for JT9D engines
                        f_JT9D = f_JT9D * 1.04; %ratio to adapt
                        f_JT9D = f_JT9D + addF;
                        T_SLS = 45550;
                        T_m = @(x) (37.707*x^2 - 47.762*x + 45.458)*1000;
                    end

                    WS_TO = WS_L/(1-X*f_JT9D); %wing loading takeoff
                    WS_IC = WS_TO * 0.965; %wing loading initial cruise
                    Cl_IC = WS_IC/(1481*PR*Mcruise^2); %cl initial cruise
                    err = abs(Cl-Cl_IC);
                    if Cl_IC > Cl
                        Cl = Cl + 0.0001;
                    else
                        Cl = Cl - 0.0001;
                    end
                end

                %TOFL W/T CALC

                Kw = syms1(TOFL/1000);
                WT_Vlo = Kw*sigma*Cl_TO/WS_TO; %thrust loading at .7Vlo
                Vlo = 1.2*sqrt((296*WS_TO)/(sigma*Cl_TO)); %1.2Vstall
                Mlo = Vlo/661/sqrt(sigma); %Mach at Vlo
                Mlo7 = Mlo*0.7; %Mach at .7Vlo

                WT = WT_Vlo * (T_m(Mlo7)/T_SLS) + addW; %thrust loading

                %WEIGHT CALC
                taperRatio = 0.35;
                Kw = (wingEngines+(fuselageEngines*1.03))/engineCount;
                n = 1.5*2.5;
                tcBar = tc+0.03;
                W_wing = (0.00945*(AR^0.8)*((1+taperRatio)^0.25)*(Kw)*(n^.5)); %numerator W_wing
                W_wing = W_wing/((tcBar^0.4)*cosd(sweep)*(WS_TO^0.695));
                %fuselage

                l = 33.2 + 3.76*PAX/abreast; %fuselage length
                d = 1 + 1.75*abreast + 1.58*aisle; %fuselage diameter
                Kf = 11.5;
                W_fuse = 0.6727*Kf*(l^.6)*(d^.72)*(n^.3);
                %landing gear
                W_lg = 0.04;
                %nacelle and pylon
                W_np = 0.0555/(WT);
                %tail surface
                Kts = (0.17*wingEngines + 0.25*fuselageEngines)/engineCount;
                W_ts = Kts;
                W_wingTS = (W_ts + 1)*W_wing;
                %power plant
                W_pp = 1/(3.58*WT);
                if engineAdv == true
                    W_pp = W_pp*1.1;
                end
                %fuel
                W_f = 1.0275 * f_JT9D;
                %payload
                W_pl = 215*PAX + W_cargo;
                %fixed equip
                Nfc = 2; %flight crew
                Nstew = PAX/50;
                W_fe = 132*PAX + 300*engineCount + 260*Nfc + 170*Nstew;
                %solve
                if composite == true
                    W_wingTS = W_wingTS * 0.7;
                    W_fuse = W_fuse *0.85;
                    W_fe = 132*PAX + 0.9*300*engineCount + 260*Nfc + 170*Nstew;
                    W_np = W_np * 0.8;
                elseif hybrid == true
                    W_wingTS = W_wingTS * 0.7;
                    W_np = W_np * 0.8;
                end
                wA = W_wingTS;
                wB = W_fuse;
                wC = W_lg + W_np + W_pp + W_f + 0.035;
                wC = wC - 1;
                wD = W_pl + W_fe;
                % W_TO = [25000:1000:600000];
                % y = wA*W_TO.^1.195 + wB*W_TO.^0.235 + wC*W_TO + wD;
                % plot(W_TO,y)

                f1 = @(WTO) wA*WTO.^1.195 + wB*WTO.^0.235 + wC*WTO + wD;
                WT_O = W_cargo*20; errW = 5000;
                while(errW>300)
                    WT_O = WT_O + errW/2;
                    errW = f1(WT_O);
                end
                %sizing
                S = WT_O / WS_TO; %wing area
                b = sqrt(AR * S); %span
                ulclc = S/b;
                T = WT_O/WT;
                Te = T/engineCount;

                %DRAG
                rnk = 0.5*sqrt(1.4*1718*394.1)/4.06E-4; % reynolds num per ft
                %Cf
                Cf_wing = RNestimateCf(rnk*ulclc)*10^-3; %skin friction coeff for wing, rnk * avg chord length
                Cf_fuse = RNestimateCf(rnk*l)*10^-3; %skin friction for fuselage, rnk * fuselage length
                %wing
                Swet_wing = 2*1.02*(S-d*30);
                Z = (2-0.5^2)*cosd(sweep)/sqrt(1-0.5^2*cosd(sweep)^2); %shevell plots
                Kw = 1 + Z*tc + 100*tc^4; %shevell plots
                f_wing = Kw*Cf_wing*Swet_wing; %flat plate drag area of wing
                %fuselage
                Swet_fuse = 0.9*pi*d*l;
                Kf = bodyForm(l/d);
                f_fuse = Kf*Cf_fuse*Swet_fuse;
                %tail
                f_ts = 0.38*f_wing;
                %nacelles
                Knac = 1.25;
                Swet_nacelles = 2.1*sqrt(Te)*engineCount;
                f_nac = Knac*Cf_wing*Swet_nacelles;
                %pylons
                f_pyl = 0.2*f_nac;
                %total
                f_total = f_wing + f_fuse + f_ts +f_nac + f_pyl;
                f_total = f_total * 1.06; %total flat plate drag
                C_do = f_total/S; %lift independant drag or form drag
                e = 1/(1.035+0.38*C_do*pi*AR); %oswald efficiency factor

                %CLIMB
                W_avgClimb = WT_O*(1+0.965)/2; %average weight during climb
                h_avgClimb = alt*20/35;
                V_ldmax = (12.9/((f_total*e)^(1/4)))*sqrt((W_avgClimb/(0.5702*b)));
                V_climb = 1.3*V_ldmax;
                Mclimb = V_climb/614.6;

                T_reqClimb = ((0.5702*f_total*V_climb^2)/296) + (94.1/(0.5702*e))*(1/(V_climb^2))*(W_avgClimb/b)^2; %thrust required for climb
                Ta_JT9D_15kClimb = @(x) -9.40672E+00*x^3 + 2.28008E+01*x^2 - 2.48184E+01*x + 2.70604E+01; %thrust chart
                Ta_JT9D_25kClimb = @(x) -6.29451E+00*x^3 + 1.42165E+01*x^2 - 1.15655E+01*x + 1.76495E+01;
                Ta_JT8D_20kClimb = @(x) -1.7117*x^3 + 4.975*x^2 - 4.2959*x + 7.4793;

                Ta_JT9D_Climb = (Ta_JT9D_25kClimb(Mclimb) +Ta_JT9D_15kClimb(Mclimb)/2)*1000; %thrust averaged
                Ta_JT8D_Climb =  Ta_JT8D_20kClimb(Mclimb)*1000;


                if engineType == 'JT9D'
                    Ta_climb = (Te/T_SLS)*Ta_JT9D_Climb; %thrust available per engine
                    SFC15k = @(x) 0.1006*x^2 + 0.3293*x + 0.3363;
                    SFC25k = @(x) -0.0107*x^2 + 0.3823*x + 0.3396;
                    SFC_c = (SFC15k(Mclimb)+SFC25k(Mclimb))/2;
                elseif engineType == "JT8D"
                    Ta_climb = (Te/T_SLS)*Ta_JT8D_Climb;
                    SFC20k = @(x) -0.0642*x^2 + 0.354*x + 0.5775;
                    SFC_c = SFC20k(Mclimb);
                end

                Ta_climbTotal = engineCount*Ta_climb;  %total thrust available
                RoC = (101*(Ta_climbTotal-T_reqClimb)*V_climb)/W_avgClimb; %rate of climb
                TimeClimb = alt/RoC;    %time to reach cruising alt
                R_climb = V_climb*TimeClimb/60; %range covered during climb
                Wf_climb = Ta_climbTotal*SFC_c*TimeClimb/60; %fuel used during climb

                %RANGE
                Wo = WT_O - Wf_climb;
                Wi = (1-W_f)*WT_O;
                Cl_avg = ((Wo+Wi)/(2*S))/(1481*PR*Mcruise^2);
                C_di = (Cl_avg^2)/(pi*AR*e); %induced drag or drag due to lift
                delta_Cdc = 0.001;
                C_d = C_di + C_do + delta_Cdc; %total drag coeff
                L_d = Cl_avg/C_d; %L/d ratio
                T_reqCruise = ((Wo+Wi)/2)/L_d;  %thrust required during cruise

                if engineType == 'JT9D'
                    T_reqJT9D = (T_reqCruise*T_SLS/Te)/engineCount; %thrust required per JT9D engine
                    SFC_35k = @(x) -0.0022*x^3 + 0.0444*x^2 - 0.3063*x + 1.3202; %@ 35k ft and M =0.82;
                    SFC_Cruise = SFC_35k(T_reqJT9D/1000);
                    if engineAdv 
                        SFC_Cruise = SFC_Cruise * 0.9;
                    end
                elseif engineType == 'JT8D'
                    T_reqJT8D = (T_reqCruise*T_SLS/Te)/engineCount;
                    SFC_35k = @(x) -0.1117*x^3 + 0.5506*x^2 - 0.3609*x + 0.8716;
                    SFC_Cruise = SFC_35k(T_reqJT8D/1000);
                end

                R_cruise = (Vcruise/SFC_Cruise)*L_d*log(Wo/Wi);
                R = R_climb + R_cruise; %actual all out range
                W_fuel = WT_O - Wi;
                errRange = abs(Rao - R);
                if Rao - R>0
                    addF = addF + 0.0001;
                    err = 1; %re-run cl loop
                elseif R - Rao >0
                    addF = addF - 0.0001;
                    err = 1; %re-run cl loop
                end
            end

            %THRUST ON TOP OF CLIMB
            Cl_top = (Wo/S)/(1481*0.236*Mcruise^2);
            Cdi_top = (Cl_top^2)/(pi*AR*e);
            Cd_top = C_do + Cdi_top + delta_Cdc;
            LD_top = Cl_top/Cd_top;
            Treq_top = (Wo/LD_top)/engineCount; %thrust required per engine
            Treq_eng = Treq_top * T_SLS/Te;
            Tavail_35k = 10000; %thrust avail at 35k
            errorWeight = Tavail_35k - Treq_eng;
            if errorWeight<0  %check
                addW = addW - 1;
                errRange = 1000;
            end
        end

        %Direct operating costs%
        CAB = range*1.15;%scheduled distance no reserves, statute miles
        Tgm = 0.25; %ground maneuver time in hours + takeoff (one min)
        Tcl = 0.18; %time for climb
        Tam = 0.1; %time for air maneuvers
        Tcr = ((CAB*1.02+20)-(R_climb*1.15))/(Vcruise*1.15); %time in cruise
        Vb = CAB/(Tgm+Tcl+Tam+Tcr); %block speed in mph
        Tb = (Tgm+Tcl+Tam+Tcr); %block time
        Fb = Wf_climb + T_reqCruise*SFC_Cruise*(Tcr+.1); %block fuel

        %FLIGHT COSTS
        payload = ((205*PAX + PAX*50)+W_cargo)/2000;
        dpbh = 40.83 + 17.849*(Vcruise*1.15*(WT_O/(10^5)))^0.3;
        Ctm_crew = dpbh/(payload*Vb);
        Ctm_fuel = (1.02*Fb*.0438 + engineCount*2.15*Tb*.135)/(CAB*payload);

        %INSURANCE COST AND HULL COST
        Wa = WT_O - W_fuel - payload - W_pp*WT_O;%airframe weight
        Ca = 87.5*Wa + 2.4*10^6; %cost of airframe
        Ce = 590000 + 16*Te; %cost of each engine
        if engineAdv == true
            Ce = Ce*1.1;
        end
        Ct = Ca + engineCount*Ce; %total cost of airplane
        U = 630 + 4000/(1+(1/(Tb+.5))); %utilization
        Ctm_ins = 0.01*Ct/(U*Vb*payload);

        %MAINTANENCE COSTS
        Kfha = 4.9169*log10(Wa/1000)-6.425;
        Kfca = 0.21256*log10(Wa/1000)^3.7375;
        Tf = Tb - Tgm;
        Ctm_direct = 8.6*(Kfha*Tf+Kfca)/(Vb*Tb*payload); %direct maint

        Cfha = 1.5994*(Ca/10^6) + 3.4263;
        Cfca = 1.9229*(Ca/10^6) + 2.2504;
        Ctm_material = (Cfha*Tf+Cfca)/(Vb*Tb*payload); %airframe material

        Kfhe = engineCount*(Te/1000)/(0.82715*(Te/1000)+13.639);
        Kfce = 0.2*engineCount;
        Ctm_engineLabor = 8.6*(Kfhe*Tf+Kfce)/(Vb*Tb*payload); %engine labor

        Cfhe = (28.2353*(Ce/10^6) - 6.5176)*engineCount;
        Cfce = (3.6698*(Ce/10^6) - 1.3685)*engineCount;
        Ctm_engineMat = (Cfhe*Tf+Cfce)/(Vb*Tb*payload); %engine material

        Ctm_maint = 2*(Ctm_direct+Ctm_material+Ctm_engineLabor+Ctm_engineMat); %total maint

        %depreciation 14 years to 10% value
        Ctm_dep = 1/(Vb*payload) * (Ct+0.06*(Ct-engineCount*Ce)+0.3*engineCount*Ce)/(14*U);

        %TOTAL DOC
        DOC = Ctm_fuel + Ctm_crew + Ctm_ins + Ctm_maint + Ctm_dep; %TOTAL DOC $/ton mile
        PAXDOC = DOC * payload/PAX; %DOC $/PAX mile
        DOCarray(j) =  PAXDOC;
        ARarr(j) = abreast;
        abreast = abreast + 1;
    end
    hold on
    figure(1)
    plot(ARarr, DOCarray)
    xlabel('Number of Seats Abreast')
    ylabel('Direct Operating Costs ($/Passenger Mi)')
    abreast = 5;
    if(i == 1)
        composite = false;
        aisle = 2;
    elseif i == 2
        composite = true;
        aisle = 1;
    elseif i == 3
        composite = true;
        aisle = 2;
    elseif i == 4
        composite = false;
        aisle = 1;
        hybrid = true;
    elseif i == 5
        composite = false;
        hybrid = true;
        aisle = 2;
    end
end

title('Direct Operating Costs with Airframe Materials, Seats Abreast and Aisles')
legend('1 Aisle - Aluminum','2 Aisle - Aluminum','1 Aisle - Composite','2 Aisle - Composite','1 Aisle- Hybrid','2 Aisle - Hybrid')
%print out
fprintf('Specifications\n')
fprintf('PAX: %i passengers   Cargo: %ilbs   Range: %inmi   Cruise Mach: %.3f   Airfoil: %s  TOFL: %ift\n',PAX,W_cargo,range,Mcruise,airfoilType,TOFL)
fprintf('AR: %i   Sweep Angle: %.1f   Wing Engines: %i   Fuselage Engines: %i   EngineType: %s   Fuel Consumed: %.2f\n', AR, sweep, wingEngines, fuselageEngines, engineType, X)
fprintf('Composites: %i Advanced Engines: %i\n\n',composite, engineAdv)
fprintf('Converged Cl: %.5f\n',Cl)
fprintf('Guess Fuel Fraction: %.5f\n', f_JT9D);
fprintf('Wing Loading - Takeoff: %.2f lbs/sqft\n', WS_TO)
fprintf('Wing Loading - Cruise: %.2f lbs/sqft\n', WS_IC)
fprintf('Thrust Loading: %.5f\n\n',WT)
fprintf('Takeoff Weight: %.2f lbs\n',WT_O);
fprintf('Wing Area: %.2f sqft\n', S)
fprintf('t/c: %.5f\n', tc)
fprintf('Span: %.3f ft\n\n', b)
fprintf('Flat Plate Drag: %.3f sqft\n',f_total);
fprintf('Form Drag Cdo: %.5f\n',C_do);
fprintf('Oswald Efficiency %.4f\n\n',e);
fprintf('Climb Speed: %3.1f knots\n',V_climb)
fprintf('Climb Thrust Required: %5.1f lbs\n',T_reqClimb)
fprintf('Climb Thrust Available: %5.1f lbs\n', Ta_climbTotal)
fprintf('SFC Climb: %.3f\n', SFC_c)
fprintf('Rate of Climb: %5.1f ft/min\n', RoC)
fprintf('Time and Range of Climb: %2.2f min and %3.1f nmi\n',TimeClimb,R_climb)
fprintf('Fuel Consumed During Climb: %4.1f lbs\n\n', Wf_climb)

fprintf('Cruising Initial and Final Weight: %6.1f lbs  | %6.1f lbs\n', Wo, Wi)
fprintf('Cl during Cruise: %.3f\n', Cl_avg)
fprintf('Induced Drag Cdi: %.3f\n', C_di)
fprintf('Combined Drag Coefficient: %.3f\n', C_d)
fprintf('SFC Cruise: %.3f\n',SFC_Cruise)
fprintf('Offset Fuel Fraction: %.5f\n',addF)
fprintf('Final Fuel Fraction: %.5f\n', W_f)
fprintf('Total Fuel Weight %5.1f lbs\n', W_fuel)
fprintf('Final All Out Range: %.1f nmi\n\n', R)

fprintf('Direct Operating Costs: %.3f $/ton mi\n',DOC)
fprintf('Direct Operating Costs: %.3f $/PAX\n', PAXDOC)