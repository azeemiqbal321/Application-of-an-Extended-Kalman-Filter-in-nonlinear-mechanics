clear all
close all
clc

samples = 26;  %Number of samples
t=[1:samples]; %Time series
R = 0; %Noise factor
measurements=[5,6,7,8,9,9,9,9,9,10,11,12,13,14,15,16,16,16,16,16,16,16,16,16,16,16];  %Noisy mass data simulation
x_est = 4;  %Initial estimate
g = 0.2;   %Scaling factor for measurement
h = 0.02;  %Scaling factor for measurement over time
dx = 1;  %Derivative of 2 , slope = 2 as dt=1
dt = 1
data = zeros(4,3)

data1 = alphaBetaFilter(samples,measurements, x_est, dx,dt, 0.302, 0.054)
data2 = alphaBetaFilter(samples,measurements, x_est, dx,dt, 0.546, 0.205)
 
%Plotting Data
plot(t,measurements,'ok','LineWidth',2)
hold on
plot(t,data1,'-b','LineWidth',2)
plot(t,data2,'-r','LineWidth',2)

legend({'Measurements',['g-0.302' char(10) 'h-0.054'],['g-0.546' char(10) 'h-0.205']})
legend('Location','southeast')
title('The g-h filter (minimizes transient errors)')
xlabel('Time (s)')
ylabel('Measurements')

