clear all
close all
clc

%Simulating Train movement
pos = 23*1000;
vel = 15;
dt=1;
measurements=pos + (vel*dt)+500*(randn(1,100));  %Noisy mass data simulation
samples = numel(measurements);  %Number of samples
t=[1:samples]; %Time series
R = 0; %Noise factor
x_est = 0;  %Initial estimate
g = 0.2;   %Scaling factor for measurement
h1 = 0.05;  %Scaling factor for measurement over time
h2 = 0.05;
h3 = 0.5;
dx = 1;  %Derivative of 2 , slope = 2 as dt=1
dt = 1  
data = zeros(4,3)

data1 = alphaBetaFilter(samples,measurements, x_est, 0,dt, g, h1)
 
%Plotting Data
plot(t,measurements)
% hold on
% plot(t,data1,'-b','LineWidth',2)

legend('Measurements','dx=0, h=0.05','dx=2, h=0.05','dx=2, h=0.5')
legend('Location','southeast')
title('The g-h filter (Varying the value of "h", the ringing effect)')
xlabel('Time (s)')
ylabel('Measurements')

