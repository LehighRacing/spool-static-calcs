FLB = 8775.38; %N 
    %Force normal to left bearing
FRB = 1633.67; %N 
    %Force normal to right bearing
    
L1 = 20; %mm 
    %Distance between diff axis and eccentric axis
L2 = 82; %mm 
    %Longest distance from diff axis to outside of eccentric plate
    
mu = .1; 
    %coefficient of friction between 7075 T6 Aluminum
TL = L1*FLB; %N­mm 
    %torque in left mount
TR = L1*FRB; %N­mm 
    %torque in right mount
FLP = 1.5*TL/L2/mu; %N 
    %force normal to left eccentric plate
FRP = 1.5*TR/L2/mu; %N 
    %force normal to right eccentric plate
FFL = TL/L2; %N 
    %frictional force left eccentric plate
FFR = TR/L2; %N 
    %frictional force right eccentric plate
FLB = FLP/2; %N 
    %force in each bolt on the left mount
FRB = FRP/2; %N 
    %force in each bolt on the right mount
    
fprintf('clamping force required per bolt (left mount) = %6.2f [N] \n', FLB)
fprintf('clamping force required per bolt (right mount) = %6.2f [N] \n', FRB)

K = .3; %Coefficient of friction from table 8­15 of shigley's
d = 10; %nominal diameter for an M8 bolt
TLB = K*FLB*d; %N­mm %torque required for each left bolt
TRB = K*FRB*d; %N­mm %torque required for each right bolt

fprintf('torque required per bolt (left mount) = %6.2f [N]\n', TLB)
fprintf('torque required per bolt (right mount) = %6.2f [N]\n', TRB)