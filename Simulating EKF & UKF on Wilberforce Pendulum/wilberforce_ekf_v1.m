clc; clear all; close all;

% Assign model parameters
global omega epsilon mass inertia
omega = 2.314;         % rad.s-1
epsilon = 9.27e-3;     % N
mass = 0.4905;         % kg
inertia = 1.39e-4;     % kg.m2
theta=2*pi;            % Initial angle
z=0;                   % Starting position

% Simulate system
xinit = [z;1;theta;0]; % initial value
h = 0.01; % time step
T = 40; %Final time value
time = 0:h:T; %Full time scale
[t,trueTrajectory]=ode45(@wilberforce_func,time,xinit);
obsNoise_z = 0.01; %Noise in position measurement
obsNoise_ang = 0.5; %Noise in angle measurement
obs_z = trueTrajectory(:,1); %First state is observable
obs_z = obs_z+obsNoise_z*randn(size(obs_z)); %Add noise
obs_ang = trueTrajectory(:,3);
obs_ang = obs_ang+obsNoise_ang*randn(size(obs_ang)); %Add noise
%plotyy(t,obs_z,t,obs_ang)

%% Linear Kalman Filter Parameters
xbar=[0;1;2*pi;0];
xbarEstimate = xbar; %Initial state
P = 0.5*eye(length(xbar)); %Initial covariance
varEstimate = diag(P); %Initial state variance
% A = [0 1 0 0;-omega^2 1 -1/2/mass*epsilon 0;0 0 0 1;-1/2/mass*epsilon 0 -omega^2 1];
% F = eye(length(xbar))+h*A; %State transition matrix
H = [1 0 0 0]; %Observation function
variance_process=5.5^2; %process variance in position z
% Q1=piecewise_white_noise(2,variance_process,h);
% Q2=piecewise_white_noise(2,variance_process,h);
Q=piecewise_white_noise(4,variance_process,h);
% Q(1:2,1:2)=Q1;
% Q(3:4,3:4)=Q2;

%Q=[h^4/4 h^3/2 0 0;h^3/2 h^2 0 0;0 0 h^4/4 h^3/2;0 0 h^3/2 h^2].*variance_process;

D = obsNoise_z;

for i = 2:length(time)
  % Defining nonlinear function
 f={@(x1,x2,x3,x4,t) (x1+x2*t);
    @(x1,x2,x3,x4,t) (x2+(-omega^2*x1-0.5/inertia*epsilon*x3)*t);
    @(x1,x2,x3,x4,t) (x3+x4*t);
    @(x1,x2,x3,x4,t) (x4+(-omega^2*x3-0.5/mass*epsilon*x1)*t)}; 
    
    % Prediction step
    xbar(1,1)=f{1}(xbar(1),xbar(2),xbar(3),xbar(4),h);
    xbar(2,1)=f{2}(xbar(1),xbar(2),xbar(3),xbar(4),h);
    xbar(3,1)=f{3}(xbar(1),xbar(2),xbar(3),xbar(4),h);
    xbar(4,1)=f{4}(xbar(1),xbar(2),xbar(3),xbar(4),h);
 A = [0 1 0 0;-omega^2 0 -0.5/inertia*epsilon 0;0 0 0 1;-0.5/mass*epsilon 0 -omega^2 0];
 F = eye(length(xbar))+h*A; %State transition matrix
 Q=piecewise_white_noise(4,variance_process,h);
 P = F*P*F' + Q;
 % Observation update
 K = (P*H')/(H*P*H'+D); % aka K = ...
 %P*H'*inv(H*P*H'+D);
 xbar = xbar + K*(obs_z(i) - H*xbar);
 P = P - K*H*P;
 xbarEstimate(:,i) = xbar;
 varEstimate(:,i) = diag(P);
 
end

figure;subplot(2,2,1);
plot(time,obs_z,'bo','markerfacecolor','b','markersize',2); hold on;
plot(time,xbarEstimate(1,:),'r','linewidth',2);
upper_z = xbarEstimate(1,:)+2*sqrt(varEstimate(1,:));
lower_z = xbarEstimate(1,:)-2*sqrt(varEstimate(1,:));
ciplot(lower_z,upper_z,time,'b',0.08);
%legend 'Noisy position' 'EKF Estimate' 'Confidence Interval'
xlabel('Time (s)')
ylabel('Position along z (m)')
% 
subplot(2,2,3);
plot(time,trueTrajectory(:,2)); hold on;
plot(time,xbarEstimate(2,:),'r','linewidth',2);
upper_z = xbarEstimate(2,:)+2*sqrt(varEstimate(2,:));
lower_z = xbarEstimate(2,:)-2*sqrt(varEstimate(2,:));
ciplot(lower_z,upper_z,time,'b',0.08);
%legend 'Simulated velocity' 'EKF Estimate' 'Confidence Interval'
xlabel('Time (s)')
ylabel('Velocity along z (m/s)')
%
subplot(2,2,2)
plot(time,obs_ang,'bo','markerfacecolor','b','markersize',2); hold on;
plot(time,xbarEstimate(3,:),'r','linewidth',3);
upper_theta = xbarEstimate(3,:)+2*sqrt(varEstimate(3,:));
lower_theta = xbarEstimate(3,:)-2*sqrt(varEstimate(3,:));
ciplot(lower_theta,upper_theta,time,'b',0.1)
%legend 'Noisy angular position' 'EKF Estimate' 'Confidence Interval'
xlabel('Time (s)')
ylabel('\theta (rad)')

subplot(2,2,4)
plot(time,trueTrajectory(:,4),'bo','markerfacecolor','b','markersize',2); hold on;
plot(time,xbarEstimate(4,:),'r','linewidth',3);
upper_theta = xbarEstimate(4,:)+2*sqrt(varEstimate(4,:));
lower_theta = xbarEstimate(4,:)-2*sqrt(varEstimate(4,:));
ciplot(lower_theta,upper_theta,time,'b',0.1)
%legend 'Simulated angular velocity' 'EKF Estimate' 'Confidence Interval'
xlabel('Time (s)')
ylabel('\omega (rad/s)')