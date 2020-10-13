clear all; close all; clc; 
data=load('data1hz.txt')
time=data(:,1);
temp=data(:,2);
plot(time,temp)
xlabel('Time (s)')
ylabel('Temperature (°C)')
title('Kalman Filter on temperature data')