clc; clear all; close all;

global gamma omega epsilon %AMP OMEG

gamma=0.3;  %Damping Coefficient
omega=-1;    %Stiffness Coefficient k1
epsilon=3;  %Stiffness Coefficient k2
% OMEG=1.25;  %Frequency of driving force
% AMP=0.5;  %Amplitude of driving force

%some parameters
Nk=10000; %Number of iterations
ndim=5;
n_meas=2; %we are only measuring the position with a displacement sensor
% % omega=sqrt(kspring/mass);
% % gamma=bspring/mass;

sensor_x_variance=0.01^2; %Position sensor var
sensor_y_variance=0.01^2; %Acceleration sensor var
dt=0.01; %time step

Q_variance=0.1; %process noise variance, assume a piecewise constant white noise
Q=piecewise_white_noise(5,Q_variance,dt);
Q=kron(Q,eye(1));

R=[sensor_x_variance 0;0 sensor_y_variance]; % Covariance for the measurement matrix

%initial state estimate
%a close to the initial measurement determines the initial state
x_initial_actual=[0; 1; 10; 10; 20];  
x=[0;1; 10; 1; 1]; %this is how the Kalman filter is initialized
P_max=10;
P=P_max*eye(5);

%choosing the sigma points, the covariance matrix is P
alpha=0.1; beta=2; kappa=3-length(x);
data=zeros(Nk,8); %placeholder for all the variables
%define the nonlinear process function called f
f={@(x1,x2,x3,x4,x5,t) (x1+x2*t);
   @(x1,x2,x3,x4,x5,t) (x2+(-x3*x2-x4^2*x1-x5*x1^3)*t)
   @(x1,x2,x3,x4,x5,t) (x3)  %gamma 
   @(x1,x2,x3,x4,x5,t) (x4)  %omega
   @(x1,x2,x3,x4,x5,t) (x5)  %epsilon
   }; 
y_sigma=zeros(ndim,2*ndim+1);

%we first simulate the behavior of the harmonic oscillator yielding true
%values
[tt,model]=ode45(@duffing_v3,dt:dt:Nk*dt,x_initial_actual(1:2));
vel = [model(1,2); model(:,2)];
tim=[tt(1,1);tt(:,1)];
difft=diff(tim);
acc = diff(vel)./difft;
acc = [acc(2);acc(2:end,1)]; 
% plot(tt,vel(2:end)); hold on;
% plot(tt,acc)

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
ZZ(1,kk)=h{1}(y_sigma(1,kk),y_sigma(2,kk),y_sigma(3,kk),y_sigma(4,kk),y_sigma(5,kk));
ZZ(2,kk)=h{2}(y_sigma(1,kk),y_sigma(2,kk),y_sigma(3,kk),y_sigma(4,kk),y_sigma(5,kk));
end

%Applying the unscented transform in the measurement space
ZZmean=sum(repmat(wm,n_meas,1).*ZZ,2); %this is our mu_z
Pz_temp=zeros(n_meas,n_meas);
for kk=1:2*ndim+1
Pz_temp=Pz_temp+wc(kk)*(ZZ(:,kk)-ZZmean)*(ZZ(:,kk)-ZZmean)';
end
Pz=Pz_temp+R;

%measurement vector z must come here
z1=model(k,1)+sqrt(sensor_x_variance)*randn; %this is the position measurement
z2=acc(k)+sqrt(sensor_y_variance)*randn; %this is the position measurement
z=[z1;z2];

y=z-ZZmean; %residual

%Finding the Kalman gain
Kg_temp=zeros(ndim,n_meas);
for kk=1:2*ndim+1
Kg_temp=Kg_temp+wc(kk)*(y_sigma(:,kk)-x)*(ZZ(:,kk)-ZZmean)';
end
Kg=Kg_temp*inv(Pz);
x=x+Kg*y;
P=P-Kg*Pz*Kg';

data(k,:)=[z(1) x(1) x(2) x(3) x(4) x(5) P(1,1) P(2,2)];
end
% % 
subplot(3,2,1)
plot(data(:,1),'ro'); hold on;
plot(data(:,2),'b-','linewidth',2);
legend('Measured displacement','Filtered displacement','location','northeast');
xlabel('t'); ylabel('x');
grid on; hold off;
% 
subplot(3,2,2)
plot(model(:,2),'ro'); hold on;
plot(data(:,3),'b-','linewidth',2);
legend('True speed','Filtered speed','location','northeast');
xlabel('t'); ylabel('speed');
grid on;

% 
subplot(3,2,3)
plot(data(:,5),'b-','linewidth',2); hold on;
%plot(ones(1,numel(data(:,4)).*omega^2),'r-','linewidth',2)
xlabel('t'); ylabel('k1');
grid on; 

subplot(3,2,4)
plot(data(:,6),'b-','linewidth',2); hold on;
%plot(ones(1,numel(data(:,4)).*epsilon),'r-','linewidth',2)
xlabel('t'); ylabel('k2');
grid on; 

subplot(3,2,[5 6])
plot(data(:,4),'b-','linewidth',2); hold on;
%plot(ones(1,numel(data(:,4)).*gamma),'r-','linewidth',2)
xlabel('t'); ylabel('gamma');
grid on; 

% % % %%%%%%%%%%
% % % 
% % % 
