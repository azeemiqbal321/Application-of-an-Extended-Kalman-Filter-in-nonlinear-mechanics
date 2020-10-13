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
xinit = [0;0;2*pi;0]; % initial value
h = 0.02; % time step
T = 100; %Final time value
time = 0:h:T; %Full time scale
[t,trueTrajectory]=ode45(@wilberforce_func,time,xinit);
obsNoise_z = 0.01; %Noise in position measurement
obsNoise_ang = 0.1; %Noise in angle measurement
obs_z = trueTrajectory(:,1); %First state is observable
obs_z = obs_z+obsNoise_z*randn(size(obs_z)); %Add noise
obs_ang = trueTrajectory(:,3);
obs_ang = obs_ang+obsNoise_ang*randn(size(obs_ang)); %Add noise
%plotyy(t,obs_z,t,obs_ang)

%% Extended Kalman Filter Parameters
xbar=[0;1;2*pi;1;1;1;1];
xbarEstimate = xbar; %Initial state
P = 0.1*eye(length(xbar)); %Initial covariance
varEstimate = diag(P); %Initial state variance
H = [1 0 0 0 0 0 0;0 0 1 0 0 0 0]; %Observation function
Q=zeros(7,7);
variance_process=2.5^2; %process variance in position z
% variance_process_theta=0.2^2; %process variance in position z
%Q=piecewise_white_noise(7,variance_process,h);
%Q2=piecewise_white_noise(2,variance_process_theta,h);
% Q(1:2,1:2)=Q1;
% Q(3:4,3:4)=Q2;
s=0.01;
G=@(t)[t^2/2;t;t^2/2;t;0.01;0.01;0.01];
Gc=@(t)[t^2/2 t t^2/2 t 0.01 0.01 0.01];
Q=variance_process*G(h)*Gc(h);
D = [obsNoise_z 0;0 obsNoise_ang];
% 
for i = 2:length(time)
 % Defining nonlinear function
 f={@(x1,x2,x3,x4,x5,x6,x7,t) (x1+x2*t);
    @(x1,x2,x3,x4,x5,x6,x7,t) (x2+(-x5*x1-x7/2/mass*x3)*t);
    @(x1,x2,x3,x4,x5,x6,x7,t) (x3+x4*t);
    @(x1,x2,x3,x4,x5,x6,x7,t) (x4+(-x6*x3-x7/2/inertia*x1)*t);
    @(x1,x2,x3,x4,x5,x6,x7,t) (x5);
    @(x1,x2,x3,x4,x5,x6,x7,t) (x6);
    @(x1,x2,x3,x4,x5,x6,x7,t) (x7)};
    % Prediction step
    xbar(1)=f{1}(xbar(1),xbar(2),xbar(3),xbar(4),xbar(5),xbar(6),xbar(7),h);
    xbar(2)=f{2}(xbar(1),xbar(2),xbar(3),xbar(4),xbar(5),xbar(6),xbar(7),h);
    xbar(3)=f{3}(xbar(1),xbar(2),xbar(3),xbar(4),xbar(5),xbar(6),xbar(7),h);
    xbar(4)=f{4}(xbar(1),xbar(2),xbar(3),xbar(4),xbar(5),xbar(6),xbar(7),h);
    xbar(5)=f{5}(xbar(1),xbar(2),xbar(3),xbar(4),xbar(5),xbar(6),xbar(7),h);
    xbar(6)=f{6}(xbar(1),xbar(2),xbar(3),xbar(4),xbar(5),xbar(6),xbar(7),h);
    xbar(7)=f{7}(xbar(1),xbar(2),xbar(3),xbar(4),xbar(5),xbar(6),xbar(7),h);
  A = [0 1 0 0 0 0 0;
      -xbar(5) 0 -0.5/mass*xbar(7) 0 -xbar(1) 0 -0.5/mass*xbar(3);
       0 0 0 1 0 0 0;
       -0.5/inertia*xbar(7) 0 -xbar(6) 0 0 -xbar(3) -0.5/inertia*xbar(1);
       0 0 0 0 0 0 0;
       0 0 0 0 0 0 0;
       0 0 0 0 0 0 0];
  F = eye(length(xbar))+h*A; %State transition matrix
  P = F*P*F' + Q;
 % Observation update
 K = (P*H')/(H*P*H'+D); % aka K = ...
 %P*H'*inv(H*P*H'+D);
 xbar = xbar + K*([obs_z(i);obs_ang(i)] - H*xbar);
 P = P - K*H*P;
 xbarEstimate(:,i) = xbar;
 varEstimate(:,i) = diag(P);
 
end

figure;subplot(3,2,1);
plot(time,obs_z,'bo','markerfacecolor','b','markersize',2); hold on;
plot(time,xbarEstimate(1,:),'r','linewidth',2);
set(gca,'fontsize',10)
legend 'Noisy position in z axis' 'EKF Estimate'
xlabel('Time (s)')
ylabel('z (m)')
% 
subplot(3,2,2)
plot(time,obs_ang,'bo','markerfacecolor','b','markersize',2); hold on;
plot(time,xbarEstimate(3,:),'r','linewidth',3);
set(gca,'fontsize',10)
legend 'Noisy angle' 'EKF Estimate'
xlabel('Time (s)')
ylabel('\theta (rad)')

subplot(3,2,3)
plot(time,xbarEstimate(5,:),'r','linewidth',3);
set(gca,'fontsize',10)
xlabel('Time (s)'); ylabel('\omega^2_z')

subplot(3,2,4)
plot(time,xbarEstimate(6,:),'r','linewidth',3);
set(gca,'fontsize',10)
xlabel('Time (s)'); ylabel('\omega^2_\theta')

subplot(3,2,[5 6])
plot(time,xbarEstimate(7,:),'r','linewidth',3);
set(gca,'fontsize',10)
xlabel('Time (s)'); ylabel('\epsilon')