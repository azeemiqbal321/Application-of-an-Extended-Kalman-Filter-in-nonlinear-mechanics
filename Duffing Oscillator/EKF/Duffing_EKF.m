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
yx=trueTrajectory(:,1); %True Position vector 
yv=trueTrajectory(:,2); %True Velocity vector

%% Applying Extended Kalman Filter
xbar=[0;1]; %Initializing the EKF
xbarEstimate=xbar;
P = diag([1,1]); %Initial covariance
varEstimate = diag(P); %Initial state variance
H = [1 0]; %Observation function
Q_variance=2.5^2; %process noise variance, assume a piecewise constant white noise
Q=piecewise_white_noise(2,Q_variance,dt);
D = obsNoise; %Measurement Noise 
var_c = 5^2; %control input noise

% Defining the nonlinear function
 f={@(x1,x2,t) (x1+x2*t);
    @(x1,x2,t) (x2+(-gamma*x2+omega^2*x1-beta*x1^3+AMP*cos(rho*t))*t)
   }; 

for i = 2:length(n_obs_pos)
%Prediction step
for kk=1:numel(xbar)
xbar(1)=f{1}(xbar(1),xbar(2),dt);
xbar(2)=f{2}(xbar(1),xbar(2),dt);
end
Aa = [0 1;-omega^2-3*beta*xbar(1)^2 -gamma];
F = eye(length(xbar))+Aa*dt+(Aa^2*dt^2)/factorial(2)+(Aa^3*dt^3)/factorial(3); %State transition matrix
P = diag([P(1,1),P(2,2)]);
M = [0 0;0 var_c];
V = [0 0;0 -AMP*sin(rho*dt)];
P = F*P*F'+V*M*V';
% Observation update
K=(P*H')/(H*P*H'+D);
y=n_obs_pos(i) - H*xbar; %Residual
xbar = xbar + K*(y);
P = P-K*H*P;
xbarEstimate(:,i) = xbar;
varEstimate(:,i) = diag(P);
residual(:,i)=y;
end

%% Position and velocity plot
subplot(2,1,1);
%plot(time,yx,'k','linewidth',2); % Plotting the true trajectory
plot(time,n_obs_pos,'o','markerfacecolor','bl','markersize',1); hold on;
lower_x1 = xbarEstimate(1,:)+2*sqrt(varEstimate(1,:));
upper_x1 = xbarEstimate(1,:)-2*sqrt(varEstimate(1,:));
ciplot(lower_x1,upper_x1,time,'b',0.05)
plot(time,xbarEstimate(1,:),'r','linewidth',2);
%set(gca,'fontsize',10)
legend({'Noisy angular position', 'Confidence Interval', 'EKF Estimate'},'FontSize',7)
xlabel('Time (s)')
hl = ylabel('$\theta~(rad)$');
set(hl, 'Interpreter', 'latex');
dim = [.14 .82 .1 .1];
str = 'a)';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none')

grid on;
% 
subplot(2,1,2)
plot(time,yv,'--','linewidth',1); % Plotting the true trajectory
%plot(time,n_obs_vel,'o','markerfacecolor','bl','markersize',1)
hold on;
lower_x2 = xbarEstimate(2,:)+2*sqrt(varEstimate(2,:));
upper_x2 = xbarEstimate(2,:)-2*sqrt(varEstimate(2,:));
ciplot(lower_x2,upper_x2,time,'b',0.05)
plot(time,xbarEstimate(2,:),'bl','linewidth',2);
legend({'True angular velocity','Confidence interval','EKF Estimate'},'FontSize',7)
xlabel('Time (s)')
hl = ylabel('$\dot{\theta}~(rad/s)$');
set(hl, 'Interpreter', 'latex');
grid on;

% %% RMSE plots
% for res = 1:numel(xbarEstimate(1,:))
% rmse_x1(res) = sqrt(sum(residual.^2)/res);
% end
% figure; 
% subplot(1,2,1)
% plot(rmse_x1,'bl','linewidth',2)
% xlabel('Number of iterations'); ylabel('RMSE')
% 
% for res2 = 1:numel(xbarEstimate(2,:))
% rmse_x2(res2) = sqrt(sum((xbarEstimate(2,:)-yv')/res2));
% end
% subplot(1,2,2)
% plot(rmse_x2,'bl','linewidth',2)
% xlabel('Number of iterations'); ylabel('RMSE')
% 
rmse_x1 = sqrt(sum((xbarEstimate(1,:)-yx').^2)/numel(time))
rmse_x2 = sqrt(sum((xbarEstimate(2,:)-yv').^2)/numel(time))
% % 
