clear;clc
%Tk= back tension at tooth k
%T0 = chain tension
%theta = sprocket minimum pressure angle 17 ? 64/N(?)
%N = number of teeth
%2*beta = sprocket tooth angle (360/N)
%k = the number of engaged teeth (angle of wrap N/360); round down to the
%nearest whole number to be safe

T0 = 7141.72; %chain tension (N)
N = 33; %number of teeth
theta = 17-64/N; %sprocket minimum pressure angle 17 ? 64/N
beta = 360/(2*N); %sprocket tooth angle 2*beta=(360/N)
k = 24; %Array for First Engaged tooth to Last Engaged Tooth

for i=1:k
Tk(i)=T0*(sind(theta)/sind(theta+2*beta))^(i-1);
end

%  uncheck_to_verify_numbers 
% disp(Tk)
for j=1:length(Tk)-1
F(j)=Tk(j)-Tk(j+1);
end
%  uncheck_to_verify_numbers 
% disp(F)
% Generate Plot of force per tooth
plot(F)
xlabel('Tooth #')
ylabel('Force on Tooth (N)')
title('Chain Force per Tooth')
grid on

axis([0 10 0 2.5e3])
fprintf('\n\n\n')

%Calculations for components of force with respect to global axis of car.
%Chain angle is -18.52 deg from positive x in the XZ Plane. From geometry
%the angle cahnge in tension direction per tooth is the same as the angle
%spacing of each tooth around the chain (360deg/N. From eximination of
%vector F, because the back tension becomes insignificant for teeth past 20
%only the first 20 teeth will be examined

delta= 360/N; %angle between teeth.
IntAng = -18.52; %Chain angle

for t=1:20
fprintf('Tooth # %6.0f\n',t)

Xcomp = F(t)*cosd(IntAng+delta*(t-1));
fprintf('Xcomp %6.2f\n',Xcomp); 
Zcomp=F(t)*sind(IntAng+delta*(t-1));
end