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
data = zeros(4,3)

data1 = alphaBetaFilter(samples,position, x_est, dx,dt, 0.2, h)
data2 = alphaBetaFilter(samples,position, x_est, dx,dt, 0.4, h)
data3 = alphaBetaFilter(samples,position, x_est, dx,dt, 0.8, h)
 

%Plotting Data
plot(t,position,'o')
hold on
plot(t,data1,'-b','LineWidth',2)
plot(t,data2,'-r','LineWidth',2)
plot(t,data3,'-g','LineWidth',2)
legend('Measurements','g-0.2','g-0.4','g-0.8')
legend('Location','southeast')
title('The g-h filter implementation (Effect of alpha (g))')
xlabel('Time (s)')
ylabel('Distance (cm)')

