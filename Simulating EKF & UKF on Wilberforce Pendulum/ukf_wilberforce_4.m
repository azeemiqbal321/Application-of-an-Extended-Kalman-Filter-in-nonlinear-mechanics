clc; clear all; close all;

%some physical properties
global omega gamma epsilon mass inertia theta z

omega = 2.314;         % m.s-1
gamma = 2.314;         % rad.s-1
epsilon = 9.27e-3;     % N
mass = 0.4905;         % kg
inertia = 1.39e-4;     % kg.m2
theta=2*pi;            % Initial angle
z=0;                   % Starting position

%some parameters
Nk=5000; %Number of iterations
ndim=7;
n_meas=2; %we are only measuring z with a motion sensor

sensor_x_variance=0.001^2;   %Position measurement error
sensor_y_variance=0.01^2;   %Angular measurement error

dt=0.01; %time step

variance_process_x=0.05^2; %process variance in position z
variance_process_y=0.5^2; %process variance in angle theta
Q1=piecewise_white_noise(2,variance_process_x,dt);
Q2=piecewise_white_noise(2,variance_process_y,dt);
Q=zeros(7,7);
Q(1:2,1:2)=Q1;
Q(3:4,3:4)=Q2;
Q=kron(Q,eye(1));

R=[sensor_x_variance 0;0 sensor_y_variance]; % Covariance for the measurement matrix

%initial state estimate
%a close to the initial measurement determines the initial state
x_initial_actual=[z;0;theta;0];  
x=[z;0;theta;0;10;10;5]; %this is how the Kalman filter is initialized
P_max=10;
P=P_max*eye(7);

%choosing the sigma points, the covariance matrix is P
alpha=0.1; beta=2; kappa=3-length(x);
data=zeros(Nk,10); %placeholder for all the variables

%define the nonlinear process function called f
f={@(x1,x2,x3,x4,x5,x6,x7,t) (x1+x2*t);
   @(x1,x2,x3,x4,x5,x6,x7,t) (x2+(-x5^2*x1-(0.5/mass*x7*x3*x1))*t);
   @(x1,x2,x3,x4,x5,x6,x7,t) (x3+x4*t);
   @(x1,x2,x3,x4,x5,x6,x7,t) (x4+(-x6^2*x3-(0.5/inertia*x7*x1*x3))*t);
   @(x1,x2,x3,x4,x5,x6,x7,t) (x5);
   @(x1,x2,x3,x4,x5,x6,x7,t) (x6);
   @(x1,x2,x3,x4,x5,x6,x7,t) (x7);
   };
y_sigma=zeros(ndim,2*ndim+1);

%we first simulate the behavior of the harmonic oscillator yielding true
%values
[tt,model]=ode45('wilberforce_func2',dt:dt:Nk*dt,x_initial_actual);

for k=1:Nk

%sigma points
    [chi,scalefactor,wm,wc]=vandermeer_sigma2(x,P,alpha,beta,kappa);
%predict step

for kk=1:2*ndim+1
y_sigma(1,kk)=f{1}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),chi(5,kk),chi(6,kk),chi(7,kk),dt);
y_sigma(2,kk)=f{2}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),chi(5,kk),chi(6,kk),chi(7,kk),dt);
y_sigma(3,kk)=f{3}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),chi(5,kk),chi(6,kk),chi(7,kk),dt);
y_sigma(4,kk)=f{4}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),chi(5,kk),chi(6,kk),chi(7,kk),dt);
y_sigma(5,kk)=f{5}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),chi(5,kk),chi(6,kk),chi(7,kk),dt);
y_sigma(6,kk)=f{6}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),chi(5,kk),chi(6,kk),chi(7,kk),dt);
y_sigma(7,kk)=f{7}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),chi(5,kk),chi(6,kk),chi(7,kk),dt);
end

%here we perform the unscented transform
x=sum(repmat(wm,ndim,1).*y_sigma,2);
P_temp=zeros(ndim,ndim);
for kk=1:2*ndim+1
P_temp=P_temp+wc(kk)*(y_sigma(:,kk)-x)*(y_sigma(:,kk)-x)';
end
P=P_temp+Q;

%define the nonlinear measurement function denoted h 
h={@(y1,y3) (y1);
   @(y1,y3) (y3);
    };

%Map the predicted prior to the measurement space
%the mapped values are stored in the variable ZZ
ZZ=zeros(n_meas,2*ndim+1);
for kk=1:2*ndim+1
ZZ(1,kk)=h{1}(y_sigma(1,kk),y_sigma(2,kk));
ZZ(2,kk)=h{2}(y_sigma(4,kk),y_sigma(3,kk));
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
zz(2)=model(k,3)+sqrt(sensor_y_variance)*randn; %this is the position measurement
y=zz'-ZZmean; %residual

%Finding the Kalman gain
Kg_temp=zeros(ndim,n_meas);
for kk=1:2*ndim+1
Kg_temp=Kg_temp+wc(kk)*(y_sigma(:,kk)-x)*(ZZ(:,kk)-ZZmean)';
end
Kg=Kg_temp*inv(Pz);
x=x+Kg*y;
P=P-Kg*Pz*Kg';

data(k,:)=[zz(1) zz(2) x(1) x(2) x(3) x(4) x(5) x(6) x(7) P(1,1)];
end
% 
figure; subplot(3,2,1);
plot(data(:,1),'ro'); hold on;
plot(data(:,3),'b-','linewidth',2);
legend('Measured','Filtered','location','northeast');
xlabel('t'); ylabel('z');
grid on; hold off;
% 
subplot(3,2,2)
plot(data(:,2),'ro');
xlabel('t'); ylabel('\theta');
grid on; hold on;
plot(data(:,5),'b-','linewidth',2);
legend('\theta','Filtered','location','northeast');
xlabel('t'); ylabel('\theta');

subplot(3,2,3)
plot(data(:,7),'b-','linewidth',2);
xlabel('t'); ylabel('\omega^2_z');
grid on; 

subplot(3,2,4)
plot(data(:,8),'b-','linewidth',2);
xlabel('t'); ylabel('\omega^2_\theta');
grid on;

subplot(3,2,[5 6])
plot(data(:,9),'b-','linewidth',2);
xlabel('t'); ylabel('\epsilon');
grid on;
% % %%%%%%%%%%
% % 
% % 
