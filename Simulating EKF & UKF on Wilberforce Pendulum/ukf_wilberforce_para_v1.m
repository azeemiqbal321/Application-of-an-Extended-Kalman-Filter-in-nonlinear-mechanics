clear all;
clc;

% Assign model parameters
global omega epsilon mass inertia theta z
omega = 2.314;         % rad.s-1
epsilon = 9.27e-3;     % N
mass = 0.4905;         % kg
inertia = 1.39e-4;     % kg.m2
theta=2*pi;            % Initial angle
z=0;                   % Starting position


%some parameters
Nk=5000; %Number of iterations
ndim=5;
n_meas=2; %we are only measuring the position with a displacement sensor

sensor_x_variance=0.001^2;
sensor_y_variance=0.01^2;
dt=0.01; %time step

Q=zeros(5,5);
Q_x_variance=0.005; %process noise variance, assume a piecewise constant white noise
Q_y_variance=0.01;
Q1=piecewise_white_noise(2,Q_x_variance,dt);
Q2=piecewise_white_noise(2,Q_y_variance,dt);
Q(1:2,1:2)=Q1;
Q(3:4,3:4)=Q2;
Q(5,5)=0;
R=[sensor_x_variance 0;0 sensor_y_variance]; % Covariance for the measurement matrix

%initial state estimate
%a close to the initial measurement determines the initial state
x_initial_actual=[z; 0; theta; 0];  
x=[0; 0; 2*pi; 0; 0.01]; %this is how the Kalman filter is initialized
P_max=10;
P=P_max*eye(5);

%choosing the sigma points, the covariance matrix is P
alpha=0.1; beta=2; kappa=3-length(x);
data=zeros(Nk,12); %placeholder for all the variables
%define the nonlinear process function called f
f={@(x1,x2,x3,x4,x5,t) (x1+x2*t);
   @(x1,x2,x3,x4,x5,t) (x2+(-omega^2*x1-0.5/mass*x5*x1^3)*t);
   @(x1,x2,x3,x4,x5,t) (x3+(-omega^2*x1-0.5/inertia*x5*x1^3)*t)
   @(x1,x2,x3,x4,x5,t) (x4);
   @(x1,x2,x3,x4,x5,t) (x5);
   };
y_sigma=zeros(ndim,2*ndim+1);

%we first simulate the behavior of the harmonic oscillator yielding true
%values
[tt,model]=ode45('wilberforce_func',dt:dt:Nk*dt,x_initial_actual(1:4));

for k=1:Nk

%sigma points
    [chi,scalefactor,wm,wc]=vandermeer_sigma2(x,P,alpha,beta,kappa);
%predict step

for kk=1:2*ndim+1
y_sigma(1,kk)=f{1}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),chi(5,kk),dt);
y_sigma(2,kk)=f{2}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),chi(5,kk),dt);
y_sigma(3,kk)=f{3}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),chi(5,kk),dt);
y_sigma(4,kk)=f{4}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),chi(5,kk),dt);
y_sigma(5,kk)=f{5}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),chi(5,kk),dt);
end

%here we perform the unscented transform
x=sum(repmat(wm,ndim,1).*y_sigma,2);
P_temp=zeros(ndim,ndim);
for kk=1:2*ndim+1
P_temp=P_temp+wc(kk)*(y_sigma(:,kk)-x)*(y_sigma(:,kk)-x)';
end
P=P_temp+Q;

%define the nonlinear measurement function denoted h 
h={@(y1,y2,y3,y4,y5) (y1);
   @(y1,y2,y3,y4,y5) (y2);
    };

%Map the predicted prior to the measurement space
%the mapped values are stored in the variable ZZ
ZZ=zeros(n_meas,2*ndim+1);
for kk=1:2*ndim+1
ZZ(1,kk)=h{1}(y_sigma(1,kk),y_sigma(2,kk),y_sigma(5,kk));
ZZ(2,kk)=h{2}(y_sigma(3,kk),y_sigma(4,kk),y_sigma(5,kk));
end

%Applying the unscented transform in the measurement space
ZZmean=sum(repmat(wm,n_meas,1).*ZZ,2); %this is our mu_z
Pz_temp=zeros(n_meas,n_meas);
for kk=1:2*ndim+1
Pz_temp=Pz_temp+wc(kk)*(ZZ(:,kk)-ZZmean)*(ZZ(:,kk)-ZZmean)';
end
Pz=Pz_temp+R;

%measurement vector z must come here
zz(1)=model(k,1)+sqrt(sensor_x_variance)*randn; %this is the position measurement
zz(2)=model(k,3)+sqrt(sensor_y_variance)*randn; %this is the angle measurement
y=zz'-ZZmean; %residual

%Finding the Kalman gain
Kg_temp=zeros(ndim,n_meas);
for kk=1:2*ndim+1
Kg_temp=Kg_temp+wc(kk)*(y_sigma(:,kk)-x)*(ZZ(:,kk)-ZZmean)';
end
Kg=Kg_temp*inv(Pz);
x=x+Kg*y;
P=P-Kg*Pz*Kg';

data(k,:)=[zz(1) zz(2) x(1) x(2) x(3) x(4) x(5) P(1,1) P(2,2) P(3,3) P(4,4) P(5,5)];
end
% 
figure; subplot(2,2,1);
plot(data(:,1),'ro'); hold on;
plot(data(:,3),'b-','linewidth',2);
legend('Measured displacement','Filtered displacement','location','northeast');
xlabel('t'); ylabel('z (m)');
grid on; hold off;
% 
subplot(2,2,2)
plot(data(:,2),'ro'); hold on;
plot(data(:,6),'b-','linewidth',2);
legend('Angular Displacement','Filtered','location','northeast');
xlabel('t'); ylabel('\theta (rad)');

subplot(2,2,[3 4])
plot(data(:,7),'b-','linewidth',2);
xlabel('t'); ylabel('epsilon');
grid on; 

% % %%%%%%%%%%
% % 
% % 
