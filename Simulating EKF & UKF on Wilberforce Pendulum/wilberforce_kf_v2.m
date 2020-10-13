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
xbar=[0;0;2*pi;0];
xbarEstimate = xbar; %Initial state
P = 1*eye(length(xbar)); %Initial covariance
varEstimate = diag(P); %Initial state variance
A = [0 1 0 0;-omega^2 1 -0.5/mass*epsilon 0;0 0 0 1; -0.5/mass*epsilon 0 -omega^2 1];
F = eye(length(xbar))+h*A; %State transition matrix
H = [1 0 0 0;0 0 1 0]; %Observation function
variance_process_x=1.5^2; %process variance in position z
variance_process_y=80^2; %process variance in angle theta
Q1=piecewise_white_noise(2,variance_process_x,h);
Q2=piecewise_white_noise(2,variance_process_y,h);
Q=zeros(4,4);
Q(1:2,1:2)=Q1;
Q(3:4,3:4)=Q2;

D = [obsNoise_z 0;0 obsNoise_ang];

% Ob = obsv(A,H);
% % If rank = ndim of state vector the system is observable
% rank(Ob)
% 
for i = 2:length(time)
 % Prediction step
 xbar = F*xbar;
 P = F*P*F' + Q;
 % Observation update
 K = (P*H')/(H*P*H'+D); % aka K = ...
 P*H'*inv(H*P*H'+D);
 xbar = xbar + K*([obs_z(i);obs_ang(i)] - H*xbar);
 P = P - K*H*P;
 xbarEstimate(:,i) = xbar;
 varEstimate(:,i) = diag(P);
 
end
figure('Renderer', 'painters', 'Position', [15 15 900 700])
subplot(2,1,1);
%plot(time,trueTrajectory(:,1),'k','linewidth',3); hold on;
plot(time,obs_z,'bo','markerfacecolor','b','markersize',2); hold on;
plot(time,xbarEstimate(1,:),'r','linewidth',2);
% upper_z = xbarEstimate(1,:)+2*sqrt(varEstimate(1,:));
% lower_z = xbarEstimate(1,:)-2*sqrt(varEstimate(1,:));
% ciplot(lower_z,upper_z,time,'b',0.08)
set(gca,'fontsize',10)
axis([0 T -0.1 0.1])
tix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(tix,'%.2f'))
dim = [.14 .82 .1 .1];
str = 'a)';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none')
legend 'Vertical Position' 'KF Estimate'
xlabel('Time (s)')
ylabel('z (m)')
% 
subplot(2,1,2)
%plot(time,trueTrajectory(:,3),'k','linewidth',2); hold on;
plot(time,obs_ang,'bo','markerfacecolor','b','markersize',2); hold on;
plot(time,xbarEstimate(3,:),'r','linewidth',3);
% upper_theta = xbarEstimate(3,:)+2*sqrt(varEstimate(3,:));
% lower_theta = xbarEstimate(3,:)-2*sqrt(varEstimate(3,:));
% ciplot(lower_theta,upper_theta,time,'b',0.1)
tix=get(gca,'ytick')';
axis([0 T -6 6])
set(gca,'yticklabel',num2str(tix,'%.2f'))
legend 'Angular Position' 'KF Estimate'
dim = [.14 .35 .1 .1];
str = 'b)';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none');
xlabel('Time (s)')
ylabel('\theta (rad)')