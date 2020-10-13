clc; clear all; close all;
%some physical properties
global kspring bspring mass
kspring=5;
bspring=3;
mass=10;

%some parameters
Nk=30000; %Number of iterations
ndim=3;
n_meas=1; %we are only measuring the position with a displacement sensor
omega=sqrt(kspring/mass);
gamma=bspring/mass;

sensor_x_variance=0.01^2;
dt=0.001; %time step

Q=zeros(3,3);
Q_variance=0.01; %process noise variance, assume a piecewise constant white noise
Q1=piecewise_white_noise(2,Q_variance,dt);
Q(1:2,1:2)=Q1;
Q(3:3,3:3)=0;
%Q=kron(Q,eye(2));
%size(Q)
R=diag([sensor_x_variance]); % Covariance for the measurement matrix

%initial state estimate
%a close to the initial measurement determines the initial state
x_initial_actual=[0; 1];  
x=[0;1;10]; %this is how the Kalman filter is initialized
P_max=10;
P=P_max*eye(3);

%choosing the sigma points, the covariance matrix is P
alpha=0.1; beta=2; kappa=3-length(x);
data=zeros(Nk,6); %placeholder for all the variables
%define the nonlinear process function called f
f={@(x1,x2,x3,t) (x1+x2*t);
   @(x1,x2,x3,t) (x2+(-gamma*x2-x3*x1)*t)
   @(x1,x2,x3,t) (x3)};
y_sigma=zeros(ndim,2*ndim+1);

%we first simulate the behavior of the harmonic oscillator yielding true
%values
[tt,model]=ode45('damped_ho',dt:dt:Nk*dt,x_initial_actual(1:2));

for k=1:Nk

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

data(k,:)=[z(1) x(1) x(2) x(3) P(1,1) P(2,2)];
end
% 
figure; subplot(2,2,1);
plot(tt,data(:,1),'ro'); hold on;
plot(tt,data(:,2),'b-','linewidth',2);
legend('Measured displacement','Filtered displacement','location','northeast');
xlabel('t'); ylabel('x');
grid on; hold off;
% 
subplot(2,2,2)
plot(tt,model(:,2),'ro'); hold on;
plot(tt,data(:,3),'b-','linewidth',2);
legend('True speed','Filtered speed','location','northeast');
xlabel('t'); ylabel('dx');

subplot(2,2,[3 4])
plot(tt,data(:,4),'b-','linewidth',2);
xlabel('t'); ylabel('omega^2');
grid on; 


% % %%%%%%%%%%
% % 
% % 
