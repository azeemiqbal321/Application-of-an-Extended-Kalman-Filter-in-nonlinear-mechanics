close all
clear all
clc

samples = 100;  %Number of samples
t=[1:samples]; %Time series
R = 100; %Noise factor
mass=2*t + R*(randn(1,samples)) + 5;  %Noisy mass data simulation
x_est = 5;  %Initial estimate
g = 0.2;   %Scaling factor for measurement
h = 0.01;  %Scaling factor for measurement over time
dx = 2;  %Derivative of 2 , slope = 2 as dt=1
dt = 1;

data = alphaBetaFilter(samples,mass, x_est, dx,dt, g, h) 


%Plotting Data
figure;plot(t,mass,'o')
hold on
plot(t,data,'-b','LineWidth',2)
legend('Measurements','Filter')
legend('Location','southeast')
title('The g-h filter implementation (Extreme Noise)')
xlabel('Time (s)')
ylabel('Mass (g)')

