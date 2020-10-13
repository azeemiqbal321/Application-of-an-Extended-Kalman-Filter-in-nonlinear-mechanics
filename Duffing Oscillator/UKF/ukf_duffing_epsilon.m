clc; clear all;
global gamma omega epsilon AMP OMEG

gamma=0.3;  %Damping Coefficient
omega=-1;    %Stiffness Coefficient k1
epsilon=1;  %Stiffness Coefficient k2
OMEG=1.25;  %Frequency of driving force
AMP=0.5;  %Amplitude of driving force

%some parameters
Nk=100; %Number of iterations
ndim=3;
n_meas=1; %we are only measuring the position with a displacement sensor
% % omega=sqrt(kspring/mass);
% % gamma=bspring/mass;

sensor_x_variance=0.5^2;
dt=2*pi/OMEG/Nk; %time step

Q_variance=5.2; %process noise variance, assume a piecewise constant white noise
Q=piecewise_white_noise(3,Q_variance,dt);
Q=kron(Q,eye(1));

R=diag([sensor_x_variance]); % Covariance for the measurement matrix

%initial state estimate
%a close to the initial measurement determines the initial state
x_initial_actual=[0; 1; 10; 10; 20];  
x=[0;1; 10]; %this is how the Kalman filter is initialized
P_max=10;
P=P_max*eye(3);

%choosing the sigma points, the covariance matrix is P
alpha=0.1; beta=2; kappa=3-length(x);
data=zeros(Nk,10); %placeholder for all the variables
%define the nonlinear process function called f
f={@(x1,x2,x3,t) (x1+x2*t);
   @(x1,x2,x3,t) (x2+(-gamma*x2+omega^2*x1-x3*x1^3+AMP*cos(OMEG))*t)
   @(x1,x2,x3,t) (x3)  %epsilon
   }; 
y_sigma=zeros(ndim,2*ndim+1);

time=0:dt:Nk;

%we first simulate the behavior of the harmonic oscillator yielding true
%values
[tt,model]=ode45(@duffing,time,x_initial_actual(1:2));

tic
for k=1:numel(model(:,1))

%sigma points
    [chi,scalefactor,wm,wc]=vandermeer_sigma2(x,P,alpha,beta,kappa);
%predict step

for kk=1:2*ndim+1
y_sigma(1,kk)=f{1}(chi(1,kk),chi(2,kk),chi(3,kk),dt);
y_sigma(2,kk)=f{2}(chi(1,kk),chi(2,kk),chi(3,kk),dt);
y_sigma(3,kk)=f{3}(chi(1,kk),chi(2,kk),chi(3,kk),dt);
end

%here we perform the unscented transform
x=sum(repmat(wm,ndim,1).*y_sigma,2);
P_temp=zeros(ndim,ndim);
for kk=1:2*ndim+1
P_temp=P_temp+wc(kk)*(y_sigma(:,kk)-x)*(y_sigma(:,kk)-x)';
end
P=P_temp+Q;

%define the nonlinear measurement function denoted h 
h={@(y1,y2,y3) (y1);
    };

%Map the predicted prior to the measurement space
%the mapped values are stored in the variable ZZ
ZZ=zeros(n_meas,2*ndim+1);
for kk=1:2*ndim+1
ZZ(1,kk)=h{1}(y_sigma(1,kk),y_sigma(2,kk),y_sigma(3,kk));
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

data(k,:)=[z(1) x(1) x(2) x(3) P(1,1) P(2,2) Kg(1) Kg(2) Kg(3) P(3,3)];
end

toc

% % 
subplot(2,2,1)
plot(time,data(:,1),'ro'); hold on;
plot(time,data(:,2),'b-','linewidth',2);
%set(gca,'fontsize',10)
legend({'Raw data', 'UKF Estimate'},'FontSize',7);
%axis([0 100 -5 5])
xlabel('Time (s)')
ylabel('$\theta$ (rad)','interpreter','latex')
grid on; hold off;
% 
subplot(2,2,2)
plot(time,model(:,2),'ro'); hold on;
plot(time,data(:,3),'b-','linewidth',2);
%set(gca,'fontsize',10)
legend({'Raw data', 'UKF Estimate'},'FontSize',7)
%axis([0 100 -5 5])
xlabel('Time (s)');
ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex')
grid on;

subplot(2,2,[3 4])
plot(time,data(:,4),'b-','linewidth',2); hold on;
xlabel('Time (s)');
ylabel('$\beta$','interpreter','latex')
axis([0 100 -10 10])
grid on; 

%% Kalman Gain plots
figure; subplot(3,1,1)
plot(time(1:200),data(1:200,7),'b-','linewidth',2); hold on;
ylabel('Gain $\theta$','interpreter','latex')
xlabel('Time (s)');
grid on; 

subplot(3,1,2)
plot(time(1:200),data(1:200,8),'b-','linewidth',2); hold on;
xlabel('Time (s)');
ylabel('Gain $\dot{\theta}$','interpreter','latex')
grid on; 

subplot(3,1,3)
plot(time(1:200),data(1:200,7),'b-','linewidth',2); hold on;
xlabel('Time (s)');
ylabel('Gain $\beta$','interpreter','latex')
grid on; 

%% Error covariance 
figure; subplot(3,1,1)
plot(time(1:200),data(1:200,5),'b-','linewidth',2); hold on;
xlabel('Time (s)');
ylabel('$P_{(1,1)}$','interpreter','latex')
grid on; 

subplot(3,1,2)
plot(time(1:200),data(1:200,6),'b-','linewidth',2); hold on;
xlabel('Time (s)');
ylabel('$P_{(2,2)}$','interpreter','latex')
grid on;

subplot(3,1,3)
plot(time(1:200),data(1:200,10),'b-','linewidth',2); hold on;
xlabel('Time (s)');
ylabel('$P_{(3,3)}$','interpreter','latex')
grid on;

rmse_1=sqrt(sum((data(:,4)-epsilon).^2)./numel(data(:,1)))


% % % %%%%%%%%%%
% % % 
% % % 
