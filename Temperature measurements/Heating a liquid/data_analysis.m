clear all; clc;

data=load('data.txt');

t=data(526:end,1)./60;
z=data(526:end,2);

% plot(t,z)
% xlabel('Time (s)');
% ylabel('Temperature (°C)');
% grid on;

Nk=numel(t);  %No. of iterations
dt=0.2;  %Time step;

R=0.0707;  %Variance of thermocouple data at room temperature
Q=0; %Processs noise covariance
xk=25;   %Initial estimate
pk=2;   %Initial error variance 
kfoutput=zeros(Nk,4);  %Predeclaration of matrix for data output

for k=1:Nk
    pk=pk+Q;    
    K = pk/(pk+R);  %Calculates Kalman Gain 
    pk = (1-K)*pk;  %Posteriori error covariance 
    y=z(k)-xk;  %Residual
    xk = xk + (K * y); %Output of filter 
    
    %Saving Data
    kfoutput(k,1)=K;
    kfoutput(k,2)=pk;
    kfoutput(k,3)=xk;
    kfoutput(k,4)=k;
end

%Plotting Data
%subplot(2,2,[1 2]) 
%plot(t,z,'-k'); hold on;
plot(t,kfoutput(:,3),'-g', 'linewidth', 2)
%tix=get(gca,'ytick')';
%set(gca,'fontsize',10)
% set(gca,'yticklabel',num2str(tix,'%.1f'))
xlabel('Time (s)')
ylabel('Temperature (°C)')
legend({'Noisy sensor data' 'Q=0.005' 'Q=1x10^{-6}' 'Q=0'}, 'location', 'southeast')
grid on;

% subplot(2,2,3)
% plot(kfoutput(1:1000,4),kfoutput(1:1000,2),'b', 'linewidth', 2)
% tix=get(gca,'ytick')';
% grid on;
% xlabel('Iteration')
% ylabel('Variance (P)')
% % 
% subplot(2,2,4)
% plot(kfoutput(1:1000,4),kfoutput(1:1000,1),'b', 'linewidth', 2)
% tix=get(gca,'ytick')';
% grid on;
% xlabel('Iteration')
% ylabel('Kalman gain (K)')

