clear all; close all; clc
newdata=load('newdata9.txt');
t=newdata(:,1);
dert=diff(t);
h=newdata(:,2); %Position
h = h';
varh=var(h);
v=diff(h)./diff(t); %Velocity
v = v';
varv=var(v);
a=diff(v)./diff(t(1:end-1)); %Acceleration
vara = var(a);
ndim=2;
dt=0.1; %sampling time
n_meas=2;
Nk=numel(h); %Number of iterations
z = h;
%%initial state x and its covariance P
x1=0.5; %Initial Height
x2=0;   %Initial Linear Velocity

x=[x1; x2]; %State Space

%choosing the sigma points, the covariance matrix is P
alpha=0.1; beta=2; kappa=3-length(x);
P = [10 0;0 10]; 

%dx=1; %actual step size to simulate the sensor measurements
%process model and its covariance Q
Q=kron(eye(2),10);

%Defining the function
f={@(x1,x2,t) (x1+x2*t);
   @(x1,x2,t) (x2)
   };

R=[5]; %Covariance of noise in measurement function

%figure; error_ellipse(P,x); title('Evolution of the state covariance P'); hold on;
%fprintf('measurement, position, speed, var(position), var(speed), covar(position,speed)\n');
data=zeros(Nk,2) %Placeholder for data points 
y_sigma=zeros(ndim,1*ndim+1); %Placeholder for sigma points

for k=1:Nk
%sigma points
 [chi,scalefactor,wm,wc]=vandermeer_sigma1(x,P,alpha,beta,kappa);
    
%predict step
for kk=1:2*ndim+1
 y_sigma(1,kk)=f{1}(chi(1,kk),chi(2,kk),dt);
 y_sigma(2,kk)=f{2}(chi(1,kk),chi(2,kk),dt);
end

x=sum(repmat(wm,ndim,1).*y_sigma,2);  %Mean of sigma points
P_temp=zeros(ndim,ndim);  

for kk=1:2*ndim+1
P_temp=P_temp+wc(kk)*(y_sigma(:,kk)-x)*(y_sigma(:,kk)-x)';
end
P=P_temp+Q;   %Covariance of sigma points

%define the nonlinear measurement function denoted h 
h={@(y1,y2) (y1);};

%Map the predicted prior to the measurement space
%the mapped values are stored in the variable ZZ
ZZ=zeros(n_meas,1*ndim+1);
for kk=1:2*ndim+1
ZZ(1,kk)=h{1}(y_sigma(1,kk),y_sigma(2,kk));
end

%Applying the unscented transform in the measurement space
ZZmean=sum(repmat(wm,n_meas,1).*ZZ,2);
Pz_temp=zeros(n_meas,n_meas);

for kk=1:2*ndim+1
Pz_temp=Pz_temp+wc(kk)*(ZZ(:,kk)-ZZmean)*(ZZ(:,kk)-ZZmean)';
end

Pz=Pz_temp+R;

%measurement vector z must come here


y=z-ZZmean; %residual

%Finding the Kalman gain
Kg_temp=zeros(ndim,n_meas);

for kk=1:2*ndim+1
Kg_temp=Kg_temp+wc(kk)*(y_sigma(:,kk)-x)*(ZZ(:,kk)-ZZmean)';
end

Kg=Kg_temp*inv(Pz);
x=x+Kg*y;
P=P-Kg*Pz*Kg';

%data(kk,:)=[z(kk) z(kk,2)];
end

figure; 
plot(t,data(:,1)); hold on;
% plot(data(:,1),data(:,2),'ro'); hold on;
% plot(data(:,3),data(:,4),'b-','MarkerSize',10,'linewidth',2);
% legend('Measured (pos)','Filtered (pos)','location','northwest');
% xlabel('x'); ylabel('y');
% grid on; hold off;
