
% ======================================================================== 
% Left Carrier Plate Forces

clear, clc, format compact, format shortg, %close all 
fprintf('============================================================== =======\n');
fprintf('Calculations for Left Spool Bracket Forces )\n');
fprintf('ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ\n\n') 

% ========================================================================= 
% Calculating the tension in the chain 
% ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ 

T_engine = 760*1.2; % [Nm] 
    % Max engine torque (w/o SF)  
d_S = 152.4/1000; % [m] 
    % drive sprocket diameter
F_T = T_engine/d_S; % [N] 
    % max Chain tension
fprintf('Max Chain Tension was approximated to be %5.1f kN\n',F_T/1000);
fprintf('(using calculations from the X49 Drivetrain Forces V2)\n\n');

BT = 154.902; % [Nm] 
    % Max braking torque (w/o SF) used for X50

% =================================================================== ====== 
% Calculations for 2D reactions in the diff brackets (in the direction of 
% chain tension). Xญdirection is towards the front of the car , which is 
% taken to be the direction of chain tension (tension has positive xญvalue 
% but negative zญvalue). 
% ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ 

% F_DL = force on bearing from Left Diff Mount 
% F_DR = force on bearing from Right Diff Mount 

% Ls_L = length between sprocket and left diff plate from CAD [m] 
% Ls_R = length between sprocket and right diff plate from CAD [m] 

Ls_L = 31.00/1000; % Length from sprocket to left diff plate [m] 
Ls_R = 142.14/1000; % Length from sprocket to right diff plate [m] 
L_lptrp = Ls_R-Ls_L; %Distance of left mount plate from right mount plate  

% Sum of Moments about Right Diff Plate: 
% F_T*(Ls_R) + F_DL*(Ls_RญLs_L) = 0 
% solving this for the force, F_DL, gives us 
F_DL = -F_T * Ls_R / L_lptrp; 

% Sum of Forces in the x direction: 
% F_T+F_DL+F_DR=0 
% solving this for the force, F_DR, gives us 
F_DR = -F_T-F_DL; 

fprintf('The force exerted on the left Diff bracket bearing is %5.1f kN\n',F_DL/1000);
fprintf('The force exerted on the right Diff bracket bearing is %5.1f kN\n',F_DR/1000);
fprintf('From the chain, the left bracket sees %3.f%% more force than the right bracket.\n\n',abs(F_DL/F_DR*100));

% NOTE: The forces calculated above for L and R diff brackets are reactions
% from the diff plates onto the bearings. For the bracket analysis, the 
% reactions from the bearings onto the diff plates must be used. Thus the 
% force below is reversed. 

F_DL = -F_DL; 

% ========================================================================= 
% Finding forces in the members of the diff brackets/supports 
% ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ 
% Step 1: Define coordinates and angles (2ญD only) from CAD model 
% The zero point is the LHS of the Rear Lower Engine mount on the YFZ (Lug)
% current file: LR-21 Spool V2.2
% ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ 
     % Z value in model, Y value in model 
LocD = [ .16353 , -0.02526]; % [m] location of diff axis 
LocU = [ -.07536 , 0.10169]; % [m] location of upper mount on engine (below intake) 
LocL = [ 0.0 , 0.0]; % [m] location of lower motor mount (rear of engine)
LocM = [ .17149, -.10599]; % [m] location of Chassis Diff Mount 
Th_Ch = -18.52; % [บ] angle of chain with ground (and therefore diff bearing reactions) 

LocB = [ 0.20595 , .01717];
LocC = [ 0.26498, 0.03283];
% Change These Below to Visualize Geometry 
% ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ 
rDz = LocM(1) - LocD(1); % [m] 
    % distance of chassis diff mount from diff axis (xญdirection) 

rDx = LocM(2) - LocD(2); % [m] 
    % distance of chassis diff mount from diff axis (zญdirection) 
rD = sqrt(rDz^2+rDx^2); % [m] 
    % distance of top diff mount from diff axis (shortest distance) 

rTx = LocB(1) - LocD(1); % [m] 
    % rTx = LocU(1) - LocD(1); % [m] 
    % distance of top diff mount from diff axis (xญdirection) 
rTz = LocB(2) - LocD(2); % [m]
    %rTz = LocU(2) - LocD(2); % [m]
    % distance of top diff mount from diff axis (zญdirection) 
rT = sqrt(rTz^2+rTx^2); % [m] 
    % distance of top diff mount from diff axis (shortest distance) 
Th_T = atand(rTz/rTx) ; % [บ] 
    % angle of top diff mount (from positive xญdirection) 
rBx = LocM(1) - LocD(1); % [m]
    % rBx = LocL(1) - LocD(1); % [m] 
    % distance of bottom diff mount from diff axis (xญdirection) 
rBz = LocM(2) - LocD(2); % [m] 
    %rBz = LocL(2) - LocD(2); % [m] 
    % distance of bottom diff mount from diff axis (zญdirection) 
rB = sqrt(rBz^2+rBx^2); % [m] 
    % distance of bottom diff mount from diff axis (shortest distance) 
Th_B = atand(rBz/rBx); % [บ]
    % angle of bottom diff mount (from positive xญdirection) 

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

LocT = LocU; % [m] 
    % location of Top mounting point on diff bracket
LocB = LocB; % [m] 
    % location of Bottom mounting point on diff bracket 

Th_BL = atand( (LocL(2) - LocB(2)) / (LocL(1) - LocB(1))); % [บ] 
    % angle of bottom turnbuckle with ground (Theta, B to L) 
Th_DT = atand( (LocC(2) - LocB(2)) / (LocC(1) - LocB(1))); % [บ] 
    % angle of bottom turnbuckle with ground (Theta, B to L) 
    
% Step 2: Sum moments about nonญ2 force reaction to find tension in % turnbuckle 
% ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ 
% Sum moment about new upper bracket mount: 
% Force xญcomponent * delta_z height + force zญcomponent * delta_x height ... =0
F_DLcomp = [F_DL*cosd(Th_Ch), 0 , F_DL*sind(Th_Ch)]; 
    % [N] components of chain force acting on diff
    
r_DLcomp = [rTx, 0, rTx]; 
    % [m] components of "r" vector from top diff pivot to diff center M = r*F ;

syms F_BL 

F_BLcomp = [F_BL*cosd(Th_B), 0, F_BL*sind(Th_B)]; % [N]
    % symbolic components of force acting on bottom diff mount
r_BLcomp = [rTx-rBx, 0, rTz-rBz]; % [m]
    % components of "r" vector from top diff pivot to bottom diff mount 
Mtemp = cross(F_DLcomp,r_DLcomp) + cross(F_BLcomp,r_BLcomp); 

% M=[r]x[F]

F_BL = double(solve(Mtemp(2),F_BL)); 

fprintf('The (bottom) Engine Mount exerts a %5.1f kN force on the diff bracket\n',F_BL/1000) 
fprintf('at an angle of %3.1f degrees (with the positive xญaxis).\n',Th_B) 
fprintf('(For the force in the member: positive ญ> tension, negative ญ> compression\n\n') 

% Step 3: Find force acting on the Top Pivot of diff (F_TP) 
% ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ 
% Sum forces acting on left diff bracket in x and z: 

F_TPx = -(F_DL*cosd(Th_Ch) + F_BL*cosd(Th_B)); 
F_TPz = -(F_DL*sind(Th_Ch) + F_BL*sind(Th_B)); 

Th_TP = atand(F_TPz/F_TPx); F_TP = sqrt(F_TPx^2 + F_TPz^2); % force acting on upper mount 
% Used to correct sign of force for upper bracket mount: 
if F_TPx <0 && F_TPz<0 
    F_TP = -F_TP; 
elseif F_TPx <0 && F_TPz>0 
        F_TP = -F_TP; 
end

fprintf('The (top) Engine Mount exerts a %5.1f kN force on the diff bracket\n',F_TP/1000)
fprintf('at an angle of %3.1f degrees (with the positive xญaxis).\n\n',Th_TP) 

% Step 4: Find forces in TU and TL. These members go from the upper diff 
% mount to the engine. Tension = positive, compression = negative 
% ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ 
% Assign angles theta_TU & theta_TL (the angle measured from the +x axis 
% of the member from the upper diffplate and theta_TU (the angle measured 
% from the +x axis of the member from the upper diffplate pivot to the 
% upper engine mount). 

Th_TU = atand((LocC(2) - LocT(2))/(LocC(1) - LocT(1))); % [บ] 
Th_TL = atand((LocM(2) - LocT(2))/(LocM(1) - LocT(1))); % [บ] 
    % Sum forces at the top pivot and solve for the unkown tensions in members  

% Eq1: Sum Fx 0 = ญF_TP*cosd(Th_TP) + F_TU*cosd(Th_TU) + F_TL*cosd(Th_TL) 
% Eq2: Sum Fy 0 = ญF_TP*sind(Th_TP) + F_TU+sind(Th_TU) + F_TL*sind(Th_TL) 

% We have two Eq's and two unknowns (F_TU and F_TL). Manipulate to solve: 
% slack JJ or email jheffmedia@gmail.com for help as to why this works 

syms F_TU F_TL 

Eq1 = -F_TP*cosd(Th_TP) + F_TU*cosd(Th_T) + F_TL*cosd(Th_B); 
Eq2 = -F_TP*sind(Th_TP) + F_TU*sind(Th_T) + F_TL*sind(Th_B); 
tempF_TL = solve(Eq1,F_TL); 
tempF_TU = double(solve(subs(Eq2,F_TL,tempF_TL),F_TU)); 
F_TL = double(subs(tempF_TL,F_TU,tempF_TU)); F_TU = tempF_TU; 
% system of equations using temp var method to solve
% symbolic solver library has issues with diverging solutions 

fprintf('Assuming pin joints for all of the engine supports, we have that following loaded members:\n') 
fprintf('ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ\n') 
fprintf('Top Diff Mount ญ> Upper Engine Mount:\n') 
fprintf('Force in member = %5.1f kN\n\n',F_TU/1000) 
fprintf('Top Diff Mount ญ> Rear Engine Mount:\n') 
fprintf('Force in member = %5.1f kN\n\n',F_TL/1000) 
fprintf('Bottom Diff Mount ญ> Rear Engine Mount:\n') 
fprintf('Force in member = %5.1f kN\n\n',F_BL/1000) 
fprintf('(For the force in the member: positive = tension, negative = compression\n\n') 

% =================================================================== ====== 
% Plot Results for 2019 Left Diff Bracket (Juicy Jay Concept)
% ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ 
clf 
figure(1) 

title('Sprocket Carriet Plate and Forces: r_T=0.060 m, \theta_T = 105^o, r_T=0.080 m, \theta_T = ญ90^o,') 
xlabel('Global xญlocation in YFZ Assembly [m]')
ylabel('Global zญlocation in YFZ Assembly [m]') 

hold on
plot([LocM(1),LocD(1),LocB(1)],[LocM(2),LocD(2),LocB(2)],'-r')
% Lines connecting diff members are shown as dashed red

plot([LocB(1),LocC(1),LocL(1),LocB(1)],[LocB(2),LocC(2),LocL(2),LocB(2)],':b') 
% Lines connecting engine members to diff are shown as dotted black 

% Plot Chassis opening (for chain) ญ can change height but thk is correct) 
ChassisOpening = [ -.05 0.075 -.05 0.305 -.065 0.305 -.065 0.075 -.05 0.075]; 
plot(ChassisOpening(:,1),ChassisOpening(:,2),'-k') 

text(.17149, -.10599,'Chassis Location') 
legend('Diff member','Engine supports','Location','Best') 
    % ADD in force on members into legend ^^
axis([-.3 .3 -0.2 0.2]) 
    % old zero point

% Color Scheme Explanation 
text(-.3,0.15,'*Forces acting on diff are black, acting on members are blue') 
AS = 40000; % Arrow Scale used below in quiver command 
rM = .020; % aesthetic radius of mount hole for plot 
rD = .065; % aesthetic radius of Diff Diameter for plot 

% Diff text and forces % ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ 
text(LocD(1)-.12,LocD(2),'F_{Chain Force Reaction at Bearing}') 
rectangle('Position',[LocD-rD/2 rD rD], 'Curvature',1, 'edgecolor','r') 

% Applied Chain Force arrow 
quiver(LocD(1), LocD(2), F_DL*cosd(Th_Ch)/AS, F_DL*sind(Th_Ch)/AS,'k') 

% Top Pivot text and forces % ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ 
text(LocT(1)+.09,LocT(2)-.05,'F_{Top Mount on Diff} (black)') 
rectangle('Position',[LocT-rM/2 rM rM], 'Curvature',1, 'edgecolor','r') 

% Top Pivot Force on Diff 
quiver(LocT(1), LocT(2), -F_TP*cosd(Th_TP)/AS, -F_TP*sind(Th_TP)/AS,'k') 

% Upper Engine Mount Force 
quiver(LocT(1), LocT(2), F_TP*cosd(Th_TP)/AS, F_TP*sind(Th_TP)/AS,'b')

% Top Diff Pivot ญ> Rear Engine Mount arrow 
quiver(LocB(1), LocB(2), -F_BL*cosd(Th_B)/AS, -F_BL*sind(Th_B)/AS,'b') 

% Turnbuckle text and forces 
% ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ 

text(LocB(1)-.055,LocB(2)-.03,'F_{TB on diff} (black)') 
text(LocB(1)+.005,LocB(2)-.02,'F_{Turnbuckle} (blue)','Color','b') 
rectangle('Position',[LocB-rM/2 rM rM], 'Curvature',1, 'edgecolor','r') 

% Lower Turnbuckle force on diff arrow 

% quiver(LocB(1), LocB(2), F_BL*cosd(Th_BL)/AS, F_BL*sind(Th_BL)/AS,'k') 
% Lower Turnbuckle member arrow %quiver(LocB(1), LocB(2), ญF_BL*cosd(Th_BL)/AS, ญF_BL*sind(Th_BL)/AS,'b') 
% Engine mount & member text and forces
% ญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญญ 
text(LocT(1)+.015,LocT(2)+.02,'F_{Upper Member}','Color','b') 

rectangle('Position',[LocT-rM/2 rM rM], 'Curvature',1, 'edgecolor','r') 
text(LocB(1)+.015,LocB(2)-.025,'F_{Lower Member}','Color','b') 
rectangle('Position',[LocB-rM/2 rM rM], 'Curvature',1, 'edgecolor','r') 
text(LocU(1)+.02,LocU(2),'Upper Engine Mount') 
rectangle('Position',[LocU-rM/2 rM rM], 'Curvature',1, 'edgecolor','k') 
text(LocL(1)+.02,LocL(2),'Rear Engine Mount') 
rectangle('Position',[LocL-rM/2 rM rM], 'Curvature',1, 'edgecolor','k') 

hold off

