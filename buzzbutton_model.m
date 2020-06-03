Nt = 3677;                 %Step Size of time
ti = 0;                   %Initial time (sec)
tf = 14.708;                  %Final time (sec)
t = linspace(ti,tf,Nt);   %Time vector (sec)

%% geometric parameters of the paperfuge system
Rr = 0.0375e-2;       % [logspace(-8,0,9),0.0375e-2];         %Radius of string (m)
Rw = 0.508e-2;               %distance between holes on the disk (m)
Rfw = 4.1275e-2;              %Radius of the full wheel (m)
Rh = 7.46125e-2;               %Radius of handle (m)
l = 22.86e-2;              %length of string (m)
w = .00762;               %thickness of wheel (m)
rho = 2710;             %density of wheel (kg/m^2)
I = (1/2)*(rho*(pi*(Rfw^2)*w))*Rfw^2;        %[logspace(-8,-2,7),9.4144e-05];   %moment of inertia of wheel (kg m^2)
Rhw = Rw + Rh;          %constant passed into ODE

%% theoretical parameters for input force
maxForce = 40;          %maximum and minimum force applied
minForce = 0.22*maxForce;%      in the force-theta phase space (N)
minThetaStart = 35;     %thetas at which force is applied
maxThetaStart = 70;     %higher maxTheta because greater torque (~60,70)%       in force-theta phase space (degrees)

inputForceList = [minThetaStart,minForce,maxThetaStart,maxForce];

%% free fitting parameters of the system (air, string)
aR = .003;              %air resistance coefficient from the wheel (kg/m^3)
gamma = 6;              %string twisting exponent
k = 1;                  %string effective spring constant (fixed at 1) (N*m)

%% run ODE with initial conditions and plot

% x1 and v1 are initial conditions that may need to be tuned so that 
%       the ODE solution does not attenuate too quickly.
%       for example, when Rfw<.01, x1=300 works well
%                    when Rfw<.006, v1 = 7000 works well
x1 = (sind(38)*l - Rhw)/Rr; %choose x1 that corresponds to theta = 38 deg
v1 = 300;

%system output "theta" is a 2000x2 array.  
%       the first column is phi (angular position) (rad)
%       the second column is phidot (angular velocity) (rad/s)
[time,phi,phiMax,phiratio,phiCrit] = solveODE_buzzbutton(t,Rr,Rhw,l,I,Rfw,w,aR,k,inputForceList,gamma,x1,v1);

%convert angular velocity from rad/sec to rpm
phi(:,2) = (phi(:,2).*60)./(2*pi);

%% plot velocity figures
figure(1);
plot(time,abs(phi(:,2)))
ylabel('Angular Velocity (RPM)')
xlabel('Time (s)')
title('Theoretical Velocity vs Experimental Velocity')

%overaly actual data on theoretical data
hold on
plot(smooth_data(:,1),smooth_data(:,2),'r');
hold off

%to compare different values of I and Rr to the max velocity

% maxV = zeros(1,length(I));
% for i = 1:length(I)
%     [time,phi,phiMax,phiratio,phiCrit] = solveODE_buzzbutton(t,Rr(10),Rhw,l,I(i),Rfw,w,aR,k,inputForceList,gamma,x1,v1);
%     maxV(i) = max(phi(:,2));
% end
% figure(2);
% semilogx(I,maxV)
% ylabel('Maximum Velocity (rad/sec')
% xlabel('Moment of Inertia (kg*m^2)')
% title('Max Velocity vs Moment of Inertia (I)')
% 
% x1R = zeros(1,length(Rr));
% for i = 1:length(Rr)
%     x1R(i) = (sind(38)*l - Rhw)/Rr(i);
%     [time,phi,phiMax,phiratio,phiCrit] = solveODE_buzzbutton(t,Rr(i),Rhw,l,I(8),Rfw,w,aR,k,inputForceList,gamma,x1R(i),v1);
%     maxV(i) = max(phi(:,2));
% end 

% figure(3);
% semilogx(Rr,maxV)
% ylabel('Maximum Velocity (rad/sec')
% xlabel('String Radius (m)')
% title('Max Velocity vs String Radius')

%% plot torque figures
theta = zeros(length(time),1);
nowF = zeros(length(time),1);
Tinput = zeros(length(time),1);
Ttwist = zeros(length(time),1);
Tdrag = zeros(length(time),1);
ratio = zeros(length(time),1);

for i = 1:length(time)
    %equation to calculate theta based on ang displacement
    theta(i) = asin((abs(phi(i,1))*Rr+Rhw)/l);
    S = sign(phi(i,1)/phi(i,2));
    if S>0 && theta(i)<minThetaStart
        nowF(i) = 0;
    elseif S>0 && theta(i)>minThetaStart
        nowF(i) = minForce;
    elseif S<0 && theta(i)>maxThetaStart
        nowF(i) = minForce;
    elseif S<0 && theta(i)<maxThetaStart
        m = (maxForce/2)/(maxThetaStart-asin(Rhw/l));
        nowF(i) = m*(theta(i)-maxThetaStart) + maxForce;
    end
    
    %use input force to calculate input torque
    if abs(phi(i,1))<phiCrit
        Tinput(i) = -sign(phi(i,1))*2*Rr*nowF(i)*(abs(phi(i,1))*Rr + Rhw)/(l^2 - (abs(phi(i,1))*Rr + Rhw)^2)^(1/2);
    else
        Tinput(i) = -sign(phi(i,1))*2*Rr*nowF(i)*tan(thetaCrit); %inputTorque fixed with theta of thetaCrit
    end
    Ttwist(i) = -sign(phi(i,1))*(1/gamma)*(phiMax-phiCrit)^(gamma+1)*((1/((phiMax-abs(phi(i,1)))^gamma))-(1/(phiMax)^gamma));
    Tdrag(i) = -sign(phi(i,2))*aR*((((4*pi)/5)*Rfw^5)+(2*pi*w*Rfw^4))*phi(i,2)^2;
    
    %to look at gamma
    ratio(i) = phi(i,1)/phiMax;
    
end

%overlay plot of different torques --> theoretical
figure(4);
plot(time,Tinput,'b')
ylabel('Torque (N*m)')
xlabel('Time (s)')
title('Torque vs Time (Theoretical Case)')

hold on
plot(time,Tdrag,'g')
plot(time,Ttwist,'r')
hold off 

%overlay plot of different torques --> experimental

%looking at gamma
figure(5);
plot(ratio,Ttwist)
ylabel('Twist Torque (N*m)')
xlabel('phi/phiMax')

hold on
plot(phiratio*ones(size(Ttwist)),Ttwist,'k')
hold off

%% plot power figures and calculate energy
% power calculations
N = 50; %number of turns
A = pi*(.03/2)^2; % in m^2
B = 0.8; % in T
Nc = 6; %number of coils
R = 130; % in milliOhms
C = 0.011*((N*A*B*Nc)^2)/R;
Ptheory = (phi(:,2).^2).*C; %have to convert angular velocity from rad/sec to rpm
Pactual = ((smooth_data(:,2)).^2).*C;

figure(6);
plot(time,Ptheory)
ylabel('Power (W)')
xlabel('Time (s)')
title('Theoretical Power vs Experimental Power')

hold on

plot(smooth_data(:,1),Pactual,'r');

hold off

%calculate electrical energy generation based on area under the curve
Etheory = trapz(time,Ptheory);
Eactual = trapz(smooth_data(:,1),Pactual);

fprintf('The theoretical energy output is %.3f Joules. The experimental energys output is %.3f Joules\n',Etheory,Eactual);

%% examine fitting parameters by comparing the theoretical velocity to the experimental
figure(7)
plot(time,abs(phi(:,2)))
ylabel('Angular Velocity (RPM)')
xlabel('Time (s)')
title('Changing aR from .003 to 3')

aR = 3;
[time,phi,phiMax,phiratio,phiCrit] = solveODE_buzzbutton(t,Rr,Rhw,l,I,Rfw,w,aR,k,inputForceList,gamma,x1,v1);
phi(:,2) = (phi(:,2).*60)./(2*pi);

hold on
plot(time,abs(phi(:,2)),'g');
plot(smooth_data(:,1),smooth_data(:,2),'r');
hold off 
