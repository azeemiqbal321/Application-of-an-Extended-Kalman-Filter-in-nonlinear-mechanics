%some parameters
Nk=50; %Number of iterations
ndim=4;
n_meas=2;
sensor_variance=0.3^2;
velocity=1;
dt=1; %time step
Q_variance=.02; %process noise variance, assume a piecewise constant white noise
R_variance=0.3^2; %mesurement noise variance
Q1=piecewise_white_noise(2,Q_variance,dt);
Q=kron(eye(2),Q1);
R=R_variance*eye(2); % Covariance for the measurement matrix
%initial state estimate
%a close to the initial measurement determines the initial state
x_initial_actual=[0; 1; 0; 1]; %[x0; vx0; y0; vy0]; 
x=[20; 5; -30; -2]; 
P_max=10;
P=P_max*eye(4);
%choosing the sigma points, the covariance matrix is P
alpha=0.1; beta=2; kappa=3-length(x);
data=zeros(Nk,8); %placeholder for all the variables
%define the nonlinear process function called f
f={@(x1,x2,x3,x4,t) (x1+x2*t);
   @(x1,x2,x3,x4,t) (x2)
   @(x1,x2,x3,x4,t) (x3+x4*t)
   @(x1,x2,x3,x4,t) (x4)
   };
y_sigma=zeros(ndim,2*ndim+1);

for k=1:Nk

%sigma points
    [chi,scalefactor,wm,wc]=vandermeer_sigma1(x,P,alpha,beta,kappa);
%predict step

for kk=1:2*ndim+1
y_sigma(1,kk)=f{1}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),dt);
y_sigma(2,kk)=f{2}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),dt);
y_sigma(3,kk)=f{3}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),dt);
y_sigma(4,kk)=f{4}(chi(1,kk),chi(2,kk),chi(3,kk),chi(4,kk),dt);
end

x=sum(repmat(wm,ndim,1).*y_sigma,2);
P_temp=zeros(ndim,ndim);
for kk=1:2*ndim+1
P_temp=P_temp+wc(kk)*(y_sigma(:,kk)-x)*(y_sigma(:,kk)-x)';
end
P=P_temp+Q;

%define the nonlinear measurement function denoted h 
h={@(y1,y2,y3,y4) (y1);
   @(y1,y2,y3,y4) (y3)
   };
%Map the predicted prior to the measurement space
%the mapped values are stored in the variable ZZ
ZZ=zeros(n_meas,2*ndim+1);
for kk=1:2*ndim+1
ZZ(1,kk)=h{1}(y_sigma(1,kk),y_sigma(2,kk),y_sigma(3,kk),y_sigma(4,kk));
ZZ(2,kk)=h{2}(y_sigma(1,kk),y_sigma(2,kk),y_sigma(3,kk),y_sigma(4,kk));
end

%Applying the unscented transform in the measurement space
ZZmean=sum(repmat(wm,n_meas,1).*ZZ,2);
Pz_temp=zeros(n_meas,n_meas);

for kk=1:2*ndim+1
Pz_temp=Pz_temp+wc(kk)*(ZZ(:,kk)-ZZmean)*(ZZ(:,kk)-ZZmean)';
end
Pz=Pz_temp+R;

%measurement vector z must come here
z=[x_initial_actual(1)+(k-1)*dt*velocity+sqrt(sensor_variance)*randn; ...
    x_initial_actual(3)+(k-1)*dt*velocity+sqrt(sensor_variance)*randn];

y=z-ZZmean; %residual

%Finding the Kalman gain
Kg_temp=zeros(ndim,n_meas);
for kk=1:2*ndim+1
Kg_temp=Kg_temp+wc(kk)*(y_sigma(:,kk)-x)*(ZZ(:,kk)-ZZmean)';
end
Kg=Kg_temp*inv(Pz);
x=x+Kg*y;
P=P-Kg*Pz*Kg';

data(k,:)=[z(1) z(2) x(1) x(3) x(2) x(4) P(1,1) P(3,3)];
end

figure; 
plot(data(:,1),data(:,2),'ro'); hold on;
plot(data(:,3),data(:,4),'b-','MarkerSize',10,'linewidth',2);
legend('Measured (pos)','Filtered (pos)','location','northwest');
xlabel('x'); ylabel('y');
grid on; hold off;

%%%%%%%%%%


