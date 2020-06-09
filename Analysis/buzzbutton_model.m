%extract real data
D = readtable('smooth_data.csv');
timeA = D{:,1};
speedA = D{:,2};
smooth_data = [timeA,speedA];

t = timeA;   %Time vector (sec)

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
        Tinput(i) = -sign(phi(i,1))*2*Rr*nowF(i)*tan(acos(2/pi)); %inputTorque fixed with theta of thetaCrit = acos(2/pi)
    end
    A = (k/gamma)*((phiMax-phiCrit)^(gamma+1));
    Ttwist(i) = -sign(phi(i,1))*A*(1/((phiMax-abs(phi(i,1)))^(gamma) - (phiMax)^-(gamma)));
    
    Tdrag(i) = -sign(phi(i,2))*aR*((((4*pi)/5)*Rfw^5)+(2*pi*w*Rfw^4))*phi(i,2)^2;
    
    %to look at gamma
    ratio(i) = phi(i,1)/phiMax;
    
end

%overlay plot of different torques --> theoretical
figure(1);
plot(time,Tinput,'b')
ylim([-0.1 0.1]);
ylabel('Torque (N*m)')
xlabel('Time (s)')
title('Torque vs Time (Theoretical Case)')

hold on
plot(time,Tdrag,'g')
plot(time,Ttwist,'r')
hold off 

legend('input Torque','drag Torque','twist Torque');

%looking at gamma
figure(2);
plot(ratio,Ttwist)
ylabel('Twist Torque (N*m)')
xlabel('phi/phiMax')

hold on
plot(phiratio*ones(size(Ttwist)),Ttwist,'k')
hold off

%% plot power figures and calculate energy
%convert angular velocity from rad/sec to rpm
phi(:,2) = (phi(:,2).*60)./(2*pi);

% power calculations
N = 50; %number of turns
A = pi*(.03/2)^2; % in m^2
B = 0.8; % in T
Nc = 6; %number of coils
R = 130; % in milliOhms
C = 0.011*((N*A*B*Nc)^2)/R;
Ptheory = (phi(:,2).^2).*C; 
Pactual = ((smooth_data(:,2)).^2).*C;

figure(3);
plot(time,Ptheory)
ylabel('Power (W)')
xlabel('Time (s)')
title('Theoretical Power vs Experimental Power')

hold on

plot(smooth_data(:,1),Pactual,'r');

hold off

%calculate energy generation based on area under the curve
Etheory = trapz(time,Ptheory);
Eactual = trapz(smooth_data(:,1),Pactual);

fprintf('The theoretical energy output is %.3f Joules. The experimental energys output is %.3f Joules\n',Etheory,Eactual);

%% plot velocity figures
figure(4);
plot(time,abs(phi(:,2)))
ylabel('Angular Velocity (RPM)')
xlabel('Time (s)')
title('Theoretical Velocity vs Experimental Velocity')

%overaly actual data on theoretical data
hold on
plot(smooth_data(:,1),smooth_data(:,2),'r');
hold off

%looking at just the 5th peak
V = readtable('smooth5.csv');
time5 = V{:,1};
speed5 = V{:,2};
smooth5 = [time5,speed5];

figure(5)
plot(smooth5(:,1),phi(1874:2315,2))
ylabel('Angular Velocity (RPM)')
xlabel('Time (s)')
title('Theoretical Velocity vs Experimental Velocity for the 5th peak')

hold on 
plot(smooth5(:,1),smooth5(:,2),'r')
hold off

legend('theoretical velocity','experimental velocity');

%to compare different parameters to the max/min velocity

% I = [logspace(-8,-2,7),9.4144e-05];
% maxV = zeros(1,length(I));
% for i = 1:length(I)
%     [time,phi,phiMax,phiratio,phiCrit] = solveODE_buzzbutton(t,Rr(10),Rhw,l,I(i),Rfw,w,aR,k,inputForceList,gamma,x1,v1);
%     maxV(i) = max(phi(:,2));
% end
% figure(6);
% semilogx(I,maxV)
% ylabel('Maximum Velocity (rad/sec')
% xlabel('Moment of Inertia (kg*m^2)')
% title('Max Velocity vs Moment of Inertia (I)')
% 

% Rr = [logspace(-8,0,9),0.0375e-2];
% x1R = zeros(1,length(Rr));
% for i = 1:length(Rr)
%     x1R(i) = (sind(38)*l - Rhw)/Rr(i);
%     [time,phi,phiMax,phiratio,phiCrit] = solveODE_buzzbutton(t,Rr(i),Rhw,l,I(8),Rfw,w,aR,k,inputForceList,gamma,x1R(i),v1);
%     maxV(i) = max(phi(:,2));
% end 

% figure(7);
% semilogx(Rr,maxV)
% ylabel('Maximum Velocity (rad/sec')
% xlabel('String Radius (m)')
% title('Max Velocity vs String Radius')

% Rfw = [logspace(-8,0,9),4.1275e-2];
% maxV = zeros(1,length(Rfw));
% minV = zeros(1,length(Rfw));
% for i = 1:length(Rfw)
%     [time,phi] = solveODE_buzzbutton(t,Rr,Rhw,l,I,Rfw(i),w,aR,k,inputForceList,gamma,x1,v1);
%     maxV(i) = max(phi(:,2));
%     minV(i) = min(phi(:,2));
% end 
% 
% figure(8);
% semilogx(Rfw,maxV)
% ylabel('Maximum Velocity (rad/sec)')
% xlabel('Radius of Disk (m)')
% title('Max Velocity vs Radius of Disk')
% 
% figure(9);
% semilogx(Rfw,minV)
% ylabel('Minimum Velocity (rad/sec)')
% xlabel('Radius of Disk (m)')
% title('Min Velocity vs Radius of Disk')
% 
% %compare distance between com and string hole (1/2 Rw) to the max velocity
% Rfw = 4.1275e-2;
% Rw = [logspace(-8,5,14),(0.508e-2)/2];
% maxV = zeros(1,length(Rw));
% for i = 1:length(Rw)
%     [time,phi] = solveODE_buzzbutton(t,Rr,(Rw(i)+Rh),l,I,Rfw,w,aR,k,inputForceList,gamma,x1,v1);
%     maxV(i) = max(phi(:,2));
% end
% figure(10);
% semilogx(Rw,maxV)
% ylabel('Maximum Velocity (rad/sec')
% xlabel('Distance Between Com and String Hole (1/2 Rw) (m)')
% title('Max Velocity vs 1/2 Rw')

%% look at phi max and phi crit experimentally and theoretically
phicrit_actual = str2double(D{(1:1:7),6});
phicrit_actual = phicrit_actual.*(2*pi); %convert from rotations to radians
thetacrit_actual = asin((phicrit_actual.*Rr)+Rhw)./l;
phimax_actual = sin(thetacrit_actual).*(l-Rhw)./Rr;
phimax_actual = phimax_actual./(2*pi); %convert back to rotations

%% examine fitting parameters by comparing the theoretical velocity to the experimental
% figure(11)
% plot(time,abs(phi(:,2)))
% ylabel('Angular Velocity (RPM)')
% xlabel('Time (s)')
% title('Changing aR from 0.003 to 3e-6')
% 
% aR = 3e-6;
% [time,phi] = solveODE_buzzbutton(t,Rr,Rhw,l,I,Rfw,w,aR,k,inputForceList,gamma,x1,v1);
% phi(:,2) = (phi(:,2).*60)./(2*pi);
% 
% hold on
% plot(time,abs(phi(:,2)),'g');
% plot(smooth_data(:,1),smooth_data(:,2),'r');
% hold off 
% legend('theoretical data, aR=0.003','theoretical data, aR=3e-6','experimental data');

