clear all
close all
clc

samples = 20;  %Number of samples
t=[1:samples]; %Time series
R = 0; %Noise factor
position=t.^2 + R*(randn(1,samples));  %Noisy mass data simulation
x_est = 5;  %Initial estimate
g = 0.2;   %Scaling factor for measurement
h = 0.02;  %Scaling factor for measurement over time
dx = 2;  %Derivative of 2 , slope = 2 as dt=1
dt = 1
data = zeros(4,4)

for k=1:samples
   %prediction step
   x_pred = x_est + (dx*dt);  
   
   %update step
   residual = position(k)-x_pred;
   dx = dx + h*(residual/dt)
   x_est = x_pred + g*residual; 
    
   %Saving Data
   %data(k,1)=t;
   data(k,2)=x_pred;
   data(k,3)=x_est;
   
end

%Plotting Data
plot(t,position,'o')
hold on
plot(t,data(:,3),'-b','LineWidth',2)
legend('Measurements','Filter')
legend('Location','southeast')
title('The g-h filter implementation (Constant Acceleration)')
xlabel('Time (s)')
ylabel('Position (cm)')

