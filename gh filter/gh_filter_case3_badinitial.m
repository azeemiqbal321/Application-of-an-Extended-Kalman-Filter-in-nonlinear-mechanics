clear all
close all
clc

samples = 100;  %Number of samples
t=[1:samples]; %Time series
R = 10; %Noise factor
mass=2*t + R*(randn(1,samples)) + 5;  %Noisy mass data simulation
x_est = 100;  %Initial estimate
g = 0.2;   %Scaling factor for measurement
h = 0.01;  %Scaling factor for measurement over time
dx = 2;  %Derivative of 2 , slope = 2 as dt=1
dt = 1
data = zeros(4,4)

for k=1:samples
   %prediction step
   x_pred = x_est + (dx*dt);  
   
   %update step
   residual = mass(k)-x_pred;
   dx = dx + h*(residual/dt)
   x_est = x_pred + g*residual; 
    
   %Saving Data
   %data(k,1)=t;
   data(k,2)=x_pred;
   data(k,3)=x_est;
   
end

%Plotting Data
plot(t,mass,'o')
hold on
plot(t,data(:,3),'-b','LineWidth',2)
legend('Measurements','Filter')
legend('Location','southeast')
title('The g-h filter implementation (Bad Initial Condition)')
xlabel('Time (s)')
ylabel('Mass (g)')

