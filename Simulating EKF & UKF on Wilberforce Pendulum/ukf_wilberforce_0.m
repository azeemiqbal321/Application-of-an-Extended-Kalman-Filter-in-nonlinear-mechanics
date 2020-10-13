clc; clear all; close all;

% Assign model parameters
global omega epsilon mass inertia theta z
omega = 2.314;         % rad.s-1
epsilon = 9.27e-3;     % N
mass = 0.4905;         % kg
inertia = 1.39e-4;     % kg.m2
theta=2*pi;            % Initial angle
z=0;                   % Starting position

xinitial = [z;0;theta;0];
dt = 0.01;     % time step
T = 30;       % Final time value
time = 0:dt:T; % Full time scale
[t,measurements]=ode45(@wilberforce_func,time,xinitial');

%some parameters
Nk=length(time); %Number of iterations
ndim=4;
n_meas=1; %we are measuring z position and angle 

sensor_x_variance=0.001^2;   %Position measurement error
variance_process_x=0.05^2; %process variance in position z
Q=piecewise_white_noise(4,variance_process_x,dt);

R=[sensor_x_variance]; % Covariance for the measurement matrix

%initial state estimate
%a close to the initial measurement determines the initial state
x_initial_actual=[z,0,theta,0];  
x=x_initial_actual'; %this is how the Kalman filter is initialized
P_max=10;
P=P_max*eye(4);

%choosing the sigma points, the covariance matrix is P
alpha=0.1; beta=2; kappa=3-length(x);
data=zeros(Nk,5); %placeholder for all the variables
%define the nonlinear process function called f
f={@(x1,x2,x3,x4,t) (x1+x2*t);
   @(x1,x2,x3,x4,t) (x2+(-omega^2*x1-0.5/mass*epsilon*x3)*t);
   @(x1,x2,x3,x4,t) (x3+x4*t)
   @(x1,x2,x3,x4,t) (x4+(-omega^2*x3-0.5/inertia*epsilon*x1)*t)};
y_sigma=zeros(ndim,2*ndim+1);

for k=1:Nk
%sigma points
    [chi,scalefactor,wm,wc]=vandermeer_sigma2(x,P,alpha,beta,kappa);

%predict step
for kk=1:2*ndim+1
y_sigma(1,kk)=f{1}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),dt);
y_sigma(2,kk)=f{2}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),dt);
y_sigma(3,kk)=f{3}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),dt);
y_sigma(4,kk)=f{4}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),dt);
end

%here we perform the unscented transform
x=sum(repmat(wm,ndim,1).*y_sigma,2);
P_temp=zeros(ndim,ndim);
for kk=1:2*ndim+1
P_temp=P_temp+wc(kk)*(y_sigma(:,kk)-x)*(y_sigma(:,kk)-x)';
end
P=P_temp+Q;

%define the nonlinear measurement function denoted h 
h={@(y1,y2,y3,y4) (y1);
    };

%Map the predicted prior to the measurement space
%the mapped values are stored in the variable ZZ
ZZ=zeros(n_meas,2*ndim+1);
for kk=1:2*ndim+1
ZZ(1,kk)=h{1}(y_sigma(1,kk),y_sigma(3,kk));
end

%Applying the unscented transform in the measurement space
ZZmean=sum(repmat(wm,n_meas,1).*ZZ,2); %this is our mu_z
Pz_temp=zeros(n_meas,n_meas);
for kk=1:2*ndim+1
Pz_temp=Pz_temp+wc(kk)*(ZZ(:,kk)-ZZmean)*(ZZ(:,kk)-ZZmean)';
end
Pz=Pz_temp+R;

%measurement vector z must come here
zz=measurements(k,1)+sqrt(sensor_x_variance)*randn; %this is the position measurement
y=zz'-ZZmean; %residual

%Finding the Kalman gain
Kg_temp=zeros(ndim,n_meas);
for kk=1:2*ndim+1
Kg_temp=Kg_temp+wc(kk)*(y_sigma(:,kk)-x)*(ZZ(:,kk)-ZZmean)';
end
Kg=Kg_temp*inv(Pz);
x=x+Kg*y;
P=P-Kg*Pz*Kg';

data(k,:)=[zz(1) x(1) x(3) P(1,1) P(3,3)];
end
% 
figure; subplot(2,1,1);
plot(t,data(:,1),'ro'); hold on;
plot(t,data(:,2),'b-','linewidth',2);
legend('Vertical displacement','Filtered','location','northeast');
xlabel('t'); ylabel('z (m)');
grid on; hold off;
tix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(tix,'%.2f'))
dim = [.14 .82 .1 .1];
str = 'a)';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none')
% 
subplot(2,1,2)
plot(t,measurements(:,3),'ro');
xlabel('t'); ylabel('\theta (rad)');
grid on; hold on;
plot(t,data(:,3),'b-','linewidth',2);
legend('Angular Displacement','Filtered','location','northeast');
dim = [.14 .35 .1 .1];
str = 'b)';
annotation('textbox',dim,'String',str,'FitBoxToText','on','LineStyle','none')
xlabel('t'); ylabel('\theta (rad)');
