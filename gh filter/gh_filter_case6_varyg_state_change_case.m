clear all
close all
clc

samples = 60;  %Number of samples
t=[1:samples]; %Time series
R = 0; %Noise factor
measurements=[5:14 14.*ones(1,50)];  %Noisy mass data simulation
x_est = 4;  %Initial estimate
g = 0.2;   %Scaling factor for measurement
h = 0.02;  %Scaling factor for measurement over time
dx = 1;  %Derivative of 2 , slope = 2 as dt=1
dt = 1
data = zeros(4,3)

data1 = alphaBetaFilter(samples,measurements, x_est, dx,dt, 0.1, h)
data2 = alphaBetaFilter(samples,measurements, x_est, dx,dt, 0.5, h)
data3 = alphaBetaFilter(samples,measurements, x_est, dx,dt, 0.9, h)
 

%Plotting Data
plot(t,measurements,'ok','LineWidth',2)
hold on
plot(t,data1,'-b','LineWidth',2)
plot(t,data2,'-r','LineWidth',2)
plot(t,data3,'-g','LineWidth',1)
legend('Measurements','g-0.1','g-0.5','g-0.9')
legend('Location','northeast')
title('The g-h filter (alpha (g) comparison in actual state change case')
xlabel('Time (s)')
ylabel('Distance (cm)')

