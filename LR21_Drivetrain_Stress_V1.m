%%Final Spool Stress
%%This analysis assumemes maximum torsion and bending. The point of of
%%analysis for Inner and Outer Diameter of the Spool Shaft is the
%%crossection at the inboard Plane of the Left (100mm id) and Right (

%Max load at left turn with full iside wheel lift
clear;clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=228*1000; % evaluated at x from Spool_Calculations2
%N*mm*****************{ENTER value at position}
%M=100*1000
T=760*1000; %N*mm
% Enter Material Yield Strength in MPa %

Sy=503; %MPa **********************{ENTER value}
Su=572; %MPa **********************{ENTER value}

%note: stress concentrations will be ignored in this analysis.
% instead, separate sections can be evaluated. Consider Expanding this in
% the future 
% Enter Inner Diameter and Outer Diameter info %
D1= 4 ; % [mm] 
    % Outer Diameter of Inner Section (inside Torsion Tube)
D2= 55 ; % [mm]
    % Outer Diameter of Outer Section (Drive Side)
d1= 24.13; % [mm]
    % Inner Diameter of Inner Section (inside Torsion Tube)
d2= 40 ; % [mm]
    % Inner Diameter of Outer Section (Drive Side)
d3 = 50;
% ******************{ENTER value}

% Calculation of Von Mises Stress

%%% from shigleys pg 326

I=(pi/64)*(D1^4-d1.^4);
    %Moment of Inertia (kg*m^2)
    
J=(pi/32)*(D1^4-d1.^4); 
    %Polar Moment of Area (kg*m^2)
    
Tau=T*(D1/2)/J; 
    %Shear Stress (MPa)
    
Sig=M*(D1/2)/I; 
    %Bending Stress (MPa)
    
Sig1=(Sig/2)+sqrt((Sig/2).^2+Tau.^2); 
    %Principle Stress 1
    
Sig2=(Sig/2)-sqrt((Sig/2).^2+Tau.^2); 
    %Principle Stress 2
    
VM=sqrt(Sig1.^2+Sig2.^2-Sig1.*Sig2); 
    %Von Mises Stress (MPa)
    
SF=Sy/VM; 
    %Static Safety Factor
    
fprintf('Static safety Factor of Concept is %6.2f\n',SF)

%%%%%%%%%%%%%%%%%%%%%%% Torsional Defection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Enter Modulus of Rigidity %%%%

G=26.9*1000; 
    %MPa***************{ENTER Value}

%%%% Length %%%%
L=170; 
    %mm******************{ENTER Value}
    
%%%%%%%%%%%%%%%%%%%%% Calculation of Deflection Angle %%%%%%%%%%%%%%%%%%%%%
theta=(T*L)/(J*G); 
    %rad
    
fprintf('Spool Angle of Twist is %6.6f (rad)\n',theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Fatigue Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sig_alt = (Sig1-Sig2)/2 ;
Sig_mean = (Sig1+Sig2)/2 ;
%%%%%%%%%%%%%%%%%%%%%%%%%Calculate Endurance Limit%%%%%%%%%%%%%%%%%%%%%%%%%

Sus=.7*Su; 
    %ferrous or not? steel­.8,other­.7
Sf=.9*Sus; 
    %torsion
Sn_prime=200; 
    %material property**************{ENTER value}

Ca=.59; 
    %for torsion from shigley 

%%% shigley (pg326) sayes to disregard if
% not purely torsional loading bc it is accounted for in distorsion
%energy theory
Ca=1; 
    %to negate its effect
Cg=.8; 
    %for 50mm<d<100mm
de=.370*D1; 
    %from shigleys: better approx then juvinall
Cg=1.51*de^(-.157);
Cs= 4.51*Su^(-.265); 
    %surface factor****************{ENTER value} from shigley table 6.2

Cs= 57.7*Su^(-.718); 
    %hot rolled

Ct=1; 
    %assumed 1

Cr=.753; 
    %at 99.9% reliability
Sn=Ca*Cg*Cs*Ct*Cr*Sn_prime;

fprintf('The Corrected Fatigue Strength is %6.3f (MPa) \n',Sn)

%%%%%%%%%%%%%%%%%%%%%%%%%% Caculate Safty Factor%%%%%%%%%%%%%%%%%%%%%%%%%%
N=(Sn*Su)/(Sig_alt*Su+Sig_mean*Sf);

fprintf('Fatigue Safety Factor is %6.2f\n',N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% New Iteration With Known Notch Parameters%%%%%%%%%%%%%%%%%%%
%The same analysis as above is conducted using Notch sensitivities. The
%above calculation give a baseline for determining the outer diameter.
%Using the Outer diameter the ratio of D/d can be determined. Also the
%fillet radius must be determined from the CAD model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%Enter Stress ConcentrationFactors%%%%%%%%%%%%%%%%%%%%
%Values from figure 4.35 in "Fundimentals of Machine Component Design" by
%(Juvinall)

Ktb=1.7; 
    %***************{ENTER Kt for bending}
Ktt=1.4; 
    %***************{ENTER Kt for torsion}

%%%%%%%%%%%%%%%%%%%%% Calculation of Von Mises Stress %%%%%%%%%%%%%%%%%%%%%
I=(pi/64)*(D1^4-d1.^4); 
    %Moment of Inertia (kg*m^2)
J=(pi/32)*(D1^4-d1.^4); 
    %Polar Moment of Area (kg*m^2)
Tau_b=Ktt*T*(D1/2)/J; 
    %Shear Stress (MPa)
Sig_b=Ktb*M*(D1/2)/I; 
    %Bending Stress (MPa)
Sig1_b=(Sig_b/2)+sqrt((Sig_b/2).^2+Tau_b.^2); 
    %Principle Stress 1
Sig2_b=(Sig_b/2)-sqrt((Sig_b/2).^2+Tau_b.^2); 
    %Principle Stress 2
VM_b=sqrt(Sig1_b.^2+Sig2_b.^2-Sig1_b.*Sig2_b); 
    %Von Mises Stress (MPa)
SF_b=Sy/VM_b; 
    %Static Safety Factor
    
fprintf('Static safety Factor w/ Stress Concentration is %6.2f\n',SF_b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Calculate Stress Concentration Factor for Fatigue%%%%%%%%%%%%%
%Esimate for q and calculation for Kf from Figure 8.24 and equation 8.2 in
%"Fundimentals of Machine Component Design" by (Juvinall)
% q= %****************{ENTER notch sensitivity factor}
% Kfb=1+(Ktb­1)*q; %for bending
% Kft=1+(Ktt­1)*q; %for torsion
%Because r>4mm q~=1 so Ktb and Ktt can be approximated by Kfb and Kft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Fatigue Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sig_alt_b=(Sig1_b-Sig2_b)/2;
Sig_mean_b=(Sig1_b+Sig2_b)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%% Caculate Safty Factor%%%%%%%%%%%%%%%%%%%%%%%%%%
N_b=(Sn*Su)/(Sig_alt_b*Su+Sig_mean_b*Sf);

nf_goodman=1/(Sig_alt_b/Sn+Sig_mean_b/Su);

nf_gerber=.5*((Su/Sig_mean_b)^2)*(Sig_alt_b/Sn)*(-1+sqrt(1+((2*Sig_mean_b*Sn)/(Su*Sig_alt_b))^2));

nf_ASME=sqrt(1/((Sig_alt_b/Sn)^2+(Sig_mean_b/Sy)^2));

fprintf('2015 Fatigue Safety Factor Equation w/ Stress Concentration is %6.2f\n',N_b)
fprintf('Goodman Safety Factor w/ Stress Concentration is %6.2f\n',nf_goodman)
fprintf('Gerber Safety Factor w/ Stress Concentration is %6.2f\n',nf_gerber)
fprintf('ASME­Elliptic Safety Factor w/ Stress Concentration is %6.2f\n',nf_ASME)