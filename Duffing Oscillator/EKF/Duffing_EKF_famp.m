clc; clear all;

%% Simulating Duffing Oscillator
global gamma omega beta AMP rho

gamma=0.3;  %Damping Coefficient
omega=-1;    %Stiffness Coefficient k1
beta=1;  %Stiffness Coefficient k2
rho=1.25;  %Frequency of driving force
AMP=0.5;  %Amplitude of driving force
T = 100; %Final time value
dt = 2*pi/rho/T; % time step
time = 0:dt:T; %Full time scale
[t,trueTrajectory]=ode45(@duffing,time,[0 1]);
obsNoise = 0.5^2; %Define observation noise level
obs_pos = trueTrajectory(:,1); %First state is observable
n_obs_pos = obs_pos+obsNoise*randn(size(obs_pos)); %Add noise
obs_vel = trueTrajectory(:,2); %Second state is observable
n_obs_vel = obs_vel+obsNoise*randn(size(obs_vel)); %Add noise
obs_epsilon=beta.*ones(numel(obs_vel),1);
yx=trueTrajectory(:,1); %Position vector 
yv=trueTrajectory(:,2); %Velocity vector
%plot(yx,yv)
%plot(n_obs,yv)

%% Applying Extended Kalman Filter
xbar=[0;1;0]; %Initializing the EKF
xbarEstimate=xbar;
P = diag([0.1,0.1,0.1]); %Initial covariance
varEstimate = diag(P); %Initial state variance
H = [1 0 0]; %Observation function
Q_variance=5.5^2; %process noise variance, assume a piecewise constant white noise
Q=piecewise_white_noise(3,Q_variance,dt);
D = obsNoise; %Measurement Noise 

for i = 2:length(n_obs_pos)
% Defining the nonlinear process function
 f={@(x1,x2,x3,t) (x1+x2*t);
    @(x1,x2,x3,t) (x2+(-gamma*x2+omega^2*x1-beta*x1^3+x3*cos(rho*t))*t)
    @(x1,x2,x3,t) (x3);
   }; 
    %Prediction step
    xbar(1)=f{1}(xbar(1),xbar(2),xbar(3),dt);
    xbar(2)=f{2}(xbar(1),xbar(2),xbar(3),dt);
    xbar(3)=f{3}(xbar(1),xbar(2),xbar(3),dt);

 Aa = [0 1 0;omega^2-3*beta*xbar(1)^2 -gamma -xbar(1)^3;0 0 0];
 F = eye(length(xbar))+Aa*dt+(Aa^2*dt^2)/factorial(2)+(Aa^3*dt^3)/factorial(3); %State transition matrix
 P = F*P*F'+Q;
 % Observation update
 K=(P*H')/(H*P*H'+D);
 xbar = xbar + K*(n_obs_pos(i) - H*xbar);
 P = P-K*H*P;
 xbarEstimate(:,i) = xbar;
 varEstimate(:,i) = diag(P);
 KalmanGain(:,i)=K;
end

%% Position and velocity plot
subplot(2,2,1);
%plot(time,yx,'k','linewidth',2); % Plotting the true trajectory
plot(time,n_obs_pos,'o','markerfacecolor','bl','markersize',1)
hold on;
plot(time,xbarEstimate(1,:),'r','linewidth',2);
set(gca,'fontsize',10)
legend({'Raw data', 'EKF Estimate'},'FontSize',7)
xlabel('Time (s)')
hl = ylabel('$\theta~(rad)$');
set(hl, 'Interpreter', 'latex');
grid on;
% 
subplot(2,2,2)
plot(time,yv,'k','linewidth',2); % Plotting the true trajectory
plot(time,n_obs_vel,'o','markerfacecolor','bl','markersize',1)
hold on;
plot(time,xbarEstimate(2,:),'bl','linewidth',2);
set(gca,'fontsize',10)
legend({'Raw data', 'EKF Estimate'},'FontSize',7)
xlabel('Time (s)')
hl = ylabel('$\dot{\theta}~(rad/s)$');
set(hl, 'Interpreter', 'latex');
grid on;
% 
subplot(2,2,[3 4])
%plot(time,epsilon.*ones(1,length(time)),'r','linewidth',2); hold on;
plot(time,xbarEstimate(3,:),'b','linewidth',2);
%legend({'True value' 'EKF Estimate'},'FontSize',7)
xlabel('Time (s)')
hl = ylabel('$\beta$');
set(hl, 'Interpreter', 'latex');
grid on;
% %%%%%%%
rmse_1 = sqrt((sum(xbarEstimate(1,:))-sum(yx)').^2/numel(yx))
rmse_2 = sqrt((sum(xbarEstimate(2,:))-sum(yv)').^2/numel(yx))
rmse_3 = sqrt((sum(xbarEstimate(3,:))-beta).^2/numel(yx))
% 
