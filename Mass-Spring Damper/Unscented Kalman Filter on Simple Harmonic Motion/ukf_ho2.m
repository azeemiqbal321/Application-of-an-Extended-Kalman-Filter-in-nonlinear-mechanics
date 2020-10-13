clc; clear all; close all;
%some physical properties
global kspring bspring mass
kspring=5;
bspring=3;
mass=10;

%some parameters
Nk=5000; %Number of iterations
ndim=4;
n_meas=1; %we are only measuring the position with a displacement sensor
omega=sqrt(kspring/mass);
gamma=bspring/mass;

sensor_x_variance=0.01^2;
dt=0.01; %time step

Q_variance=0.01; %process noise variance, assume a piecewise constant white noise
Q1=piecewise_white_noise(2,Q_variance,dt);
Q=kron(Q1,eye(2));
R=diag([sensor_x_variance]); % Covariance for the measurement matrix

%initial state estimate
%a close to the initial measurement determines the initial state
x_initial_actual=[0; 1; 10; 10];  
x=[-2;1.5; 0; 0]; %this is how the Kalman filter is initialized
P_max=10;
P=P_max*eye(4);

%choosing the sigma points, the covariance matrix is P
alpha=0.1; beta=2; kappa=3-length(x);
data=zeros(Nk,7); %placeholder for all the variables
%define the nonlinear process function called f
f={@(x1,x2,x3,x4,t) (x1+x2*t);
   @(x1,x2,x3,x4,t) (x2+(-(x3/mass)*x1-(x4/mass)*x2)*t)
   @(x1,x2,x3,x4,t) (x3)
   @(x1,x2,x3,x4,t) (x4)};
y_sigma=zeros(ndim,2*ndim+1);

time=dt:dt:Nk*dt;

%we first simulate the behavior of the harmonic oscillator yielding true
%values
[tt,model]=ode45('damped_ho',time,x_initial_actual(1:2));

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

data(k,:)=[z(1) x(1) x(2) x(3) x(4) P(1,1) P(2,2)];
end
% 
figure('Renderer', 'painters', 'Position', [15 15 900 700])
subplot(2,2,1);
plot(time,data(:,1),'ro'); hold on;
plot(time,data(:,2),'b-','linewidth',2);
tix=get(gca,'ytick')';
set(gca,'fontsize',10)
set(gca,'yticklabel',num2str(tix,'%.1f'))
axis([0 30,-1 1.5])
legend('Measured displacement','Filtered displacement','location','northeast');
xlabel('Time (s)'); ylabel('Position (m)');
grid on; hold off;
% 
subplot(2,2,2)
plot(time,model(:,2),'ro'); grid on; hold on;
plot(time,data(:,3),'b-','linewidth',2);
set(gca,'fontsize',10)
set(gca,'yticklabel',num2str(tix,'%.1f'))
axis([0 30,-1 1.5])
legend('True speed','Filtered speed','location','northeast');
xlabel('Time (s)'); ylabel('Velocity (m/s)');

subplot(2,2,3)
plot(time,data(:,4),'b-','linewidth',2);
%axis([0 30],[0 10])
set(gca,'fontsize',10)
% set(gca,'yticklabel',num2str(tix,'%.1f'))
xlabel('Time (s)'); ylabel('Spring Constant (N/m)');
grid on; 

subplot(2,2,4)
%plot(gamma.*ones(1,Nk),'r','linewidth',4); hold on;
plot(time,data(:,5),'b-','linewidth',2);
%axis([0 30],[0 10])
set(gca,'fontsize',10)
% set(gca,'yticklabel',num2str(tix,'%.1f'))
xlabel('Time (s)'); ylabel('Damping Coefficient (Ns/m)');
grid on; 
% % %%%%%%%%%%
% % 
% % 

rmse_1 = sqrt(sum((data(:,4)-kspring.*ones(numel(data(:,3)),1)).^2)./numel(data(:,1)))
rmse_2 = sqrt(sum((data(:,5)-bspring.*ones(numel(data(:,3)),1)).^2)./numel(data(:,1)))
% rmse_2 = sqrt((mean(data(:,3)-model(:,2)).^2)./numel(data(:,1)))
% rmse_3 = sqrt((mean(data(:,4)'-kspring.*ones(1,Nk)).^2)./numel(data(:,1)))
% rmse_4 = sqrt((mean(data(:,5)'-bspring.*ones(1,Nk)).^2)./numel(data(:,1)))
