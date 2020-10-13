clc; clear all; close all;

global gamma omega epsilon AMP OMEG
gamma=0.3;  %Damping Coefficient
omega=-1;    %Stiffness Coefficient k1
epsilon=1;  %Stiffness Coefficient k2
OMEG=1.25;  %Frequency of driving force
AMP=0.5;  %Amplitude of driving force

%some parameters
ndim=2;
n_meas=1; %we are only measuring the position with a displacement sensor
sensor_x_variance=0.1^2; %Measurement noise
x_initial_actual=[0;1]; 
dt = 2*pi/OMEG/100; % time step
T = 100; %Final time value
time = 0:dt:T; %Full time scale
%we first simulate the behavior of the duffing oscillator yielding true
%values
[tt,model]=ode45(@duffing,time,x_initial_actual(1:2));
Nk=numel(model(:,1)); %Number of iterations


Q_variance=0.5^2; %process noise variance, assume a piecewise constant white noise
Q=piecewise_white_noise(2,Q_variance,dt);
Q=kron(Q,eye(1));

R=diag([sensor_x_variance]); % Covariance for the measurement matrix

%initial state estimate
%a close to the initial measurement determines the initial state 
x=[0;1]; %this is how the Kalman filter is initialized
P_max=10;
P=P_max*eye(2);

%choosing the sigma points, the covariance matrix is P
alpha=0.1; beta=2; kappa=3-length(x);
data=zeros(Nk,5); %placeholder for all the variables
%define the nonlinear process function called f
f={@(x1,x2,t) (x1+x2*t);
   @(x1,x2,t) (x2+(-gamma*x2+omega^2*x1-epsilon*x1^3+AMP*cos(OMEG))*t);
   }; 
y_sigma=zeros(ndim,2*ndim+1);

for k=1:Nk
%sigma points
    [chi,scalefactor,wm,wc]=vandermeer_sigma2(x,P,alpha,beta,kappa);
%predict step
for kk=1:2*ndim+1
y_sigma(1,kk)=f{1}(chi(1,kk),chi(2,kk),dt);
y_sigma(2,kk)=f{2}(chi(1,kk),chi(2,kk),dt);
end

%here we perform the unscented transform
x=sum(repmat(wm,ndim,1).*y_sigma,2);
P_temp=zeros(ndim,ndim);
for kk=1:2*ndim+1
P_temp=P_temp+wc(kk)*(y_sigma(:,kk)-x)*(y_sigma(:,kk)-x)';
end
P=P_temp+Q;

%define the nonlinear measurement function denoted h 
h={@(y1,y2) (y1);
    };

%Map the predicted prior to the measurement space
%the mapped values are stored in the variable ZZ
ZZ=zeros(n_meas,2*ndim+1);
for kk=1:2*ndim+1
ZZ(1,kk)=h{1}(y_sigma(1,kk),y_sigma(2,kk));
end

%Applying the unscented transform in the measurement space
ZZmean=sum(repmat(wm,n_meas,1).*ZZ,2); %this is our mu_z
Pz_temp=zeros(n_meas,n_meas);
for kk=1:2*ndim+1
Pz_temp=Pz_temp+wc(kk)*(ZZ(:,kk)-ZZmean)*(ZZ(:,kk)-ZZmean)';
end
Pz=Pz_temp+R;

%measurement vector z must come here
z=model(k,1)+sqrt(sensor_x_variance)*randn; %this is the position measurement
y=z-ZZmean; %residual

%Finding the Kalman gain
Kg_temp=zeros(ndim,n_meas);
for kk=1:2*ndim+1
Kg_temp=Kg_temp+wc(kk)*(y_sigma(:,kk)-x)*(ZZ(:,kk)-ZZmean)';
end
Kg=Kg_temp*inv(Pz);
x=x+Kg*y;
P=P-Kg*Pz*Kg';

data(k,:)=[z(1) x(1) x(2) P(1,1) P(2,2)];
end
% % 
%% Position and velocity plot
figure(2);subplot(2,1,1);
% plot(tt,model(:,1),'k','linewidth',2); % Plotting the true trajectory
% hold on;
plot(tt,data(:,1),'o','markerfacecolor','bl','markersize',1); hold on;
plot(tt,data(:,2),'r','linewidth',2);
lower_x1 = data(:,2)+2*sqrt(data(:,4));
upper_x1 = data(:,2)-2*sqrt(data(:,4));
ciplot(lower_x1,upper_x1,tt,'b')
tix=get(gca,'ytick')';
set(gca,'fontsize',10)
set(gca,'yticklabel',num2str(tix,'%.1f'))
legend({'Noisy angular position' 'UKF estimate' 'Confidence interval'},'FontSize',7);
xlabel('Time (s)')
ylabel('$\theta~$(rad)','interpreter','latex')
% 
subplot(2,1,2)
plot(tt,model(:,2)); hold on;
axis([0 100 -4.0 4.0])
plot(tt,data(:,3),'b','linewidth',2);
lower_x1 = data(:,3)+2*sqrt(data(:,5));
upper_x1 = data(:,3)-2*sqrt(data(:,5));
ciplot(lower_x1,upper_x1,tt,'b')
tix=get(gca,'ytick')';
set(gca,'fontsize',10)
set(gca,'yticklabel',num2str(tix,'%.1f'))
legend({'Angular velocity' 'UKF estimate' 'Confidence interval'},'FontSize',7);
xlabel('Time (s)')
ylabel('$\dot{\theta}~$(rad/s)','interpreter','latex')

rmse_1=sqrt(sum((data(:,1)-data(:,2)).^2)./numel(data(:,1)))
rmse_2=sqrt(sum((model(:,2)-data(:,3)).^2)./numel(data(:,1)))