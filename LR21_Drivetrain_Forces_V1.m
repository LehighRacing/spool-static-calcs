%% SPOOL Calculations for Shear and Moment Diagrams 
%Using: torque from 2014 report. Using: sprocket radius and coordinates of 
%Sprocket plane, Left Bearing Plane, Right Bearing Plane (lengths can be %determined from basic coordinate math in Y) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear; 
clc 

Torque = 760; %Output Torque (N*m) 

SprocketR= 76.2/1000; %drive sprocket radius m 

S=Torque/SprocketR; %Chain Tension (N) 
fprintf('Chain Tension is %6.2f (N)\n',S)
fprintf('at angle of 11.03 deg from X axis \n\n\n')

L1=31; % [mm] 
    % From Sprocket to Center L Bearing
L2=111.4; % [mm]
    % From Center L Bearing to Center R Bearing
L3=L1+L2; % [mm]
    % From Sprocket to center R Bearing  

%%%%%%%%%%%%%%%%%Sudo­Code for Sum of Forces and Moments%%%%%%%%%%%%%%%%%%% 
% From Sum of Forces and Moments diagram (see pdf notes)
% Sum F= S ­ FL + FR = 0
% Sum M (about S)= ­FL*L1 + FR*L3 = 0 

%%%%%%%%%% Solved Equations %%%%%%%%%%
FL=S*(L3/L2); % Force in Left Bearing (N) 
FR=FL-S; % Force in Right Bearing (N)

%%%%%%%% Output %%%%%%%%
fprintf('Force in Left Bearing %6.2f (N) \n',FL) 
fprintf('Force in Right Bearing %6.2f (N) \n',FR) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%LBI=Left Bearing Inboard LBO=Left Bearing Outboard 
%RBI=Right Bearing Inboard RBO=Right Bearing Outboard 
%Sp=Sprocket 

%%%%%%%%%%%%%%%%%%%%%%% Equations for F V M diagrams%%%%%%%%%%%%%%%%%%%%%%%

l1= 23; %Sp to LBO (mm) 
l2= 18; %LBO to LBI (mm) [Left Bearing Width]
l3=95.14; %LBI to RBI (mm)
l4= 16; %RBI to RBO (mm) [Right Bearing Width]

DL = FL/l2 %L bearing distibuted load (N/mm)
DR = FR/l4 %R bearing distibuted load (N/mm) 
fprintf('\n\n') 

%V 

V1=S; %Shear (from Sp to LBO) (N) 
V2=V1-(DL*l2); %Shear (at LBI) (N) 
V3=V2+(DR*l4); %Shear (at RBO) (N)
fprintf('Shear from Sprocket to Outer Left Bearing %6.2f (kN) \n',V1/1000) 
fprintf('Shear at Inner Left Bearing %6.2f (kN) \n',V2/1000)
fprintf('Shear at Outer Right Bearing %6.2f (kN) \n',V3/1000) 
fprintf('\n\n') 

%Zero Points 
Z1=V1/DL; %Location of zero V and Max M (measured from LBO) (mm) 
Z2=l2-Z1; %Location of zero V and Max M (measured from Z1 to LBO) 
fprintf('Location of zero V and Max M (measured from LBO) %6.2f (mm) \n',Z1)
fprintf('%Location of zero V and Max M (measured from Z1 to LBO) %6.2f (mm) \n',Z2) 

%M
M1=V1*l1; %Moment at LBO (N*mm) 
M2=M1+(DL/2)*(Z1^2); %Max Moment at Z1 (N*mm) 
M3=M2-(DL/2)*(Z2^2); %Moment at LBI (N*mm) 
M4=M3+V2*l3; %Moment at RBI (N*mm) 
M5=M4-(DR/2)*(l4^2); %Moment at RBO (N*mm)

x = 23.841 ; %[mm]
    % probe distance measured from LBI 

Mx=M3-((M3-M4)/(l3))*x; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
fprintf('Moment on Outer Side Left Bearing %6.2f (N*m) \n',M1/1000)
fprintf('Maximum Moment on Left Bearing %6.2f (N*m) \n',M2/1000) 
fprintf('Moment on Inner Side Left Bearing %6.2f (N*m) \n',M3/1000) 
fprintf('Moment on Inner Side Right Bearing %6.2f (N*m) \n',M4/1000)
fprintf('Moment on Outer Side Right Bearing %6.2f (N*m) \n',M5/1000)
fprintf('Moment on Probe Distance %6.2f (N*m) \n',Mx/1000)

