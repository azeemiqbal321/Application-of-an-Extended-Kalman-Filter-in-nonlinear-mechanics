close all
clear all
clc

%Dog's initial position
t=[6:0.1:14];
x = normpdf(t,10.2,1)
plot(t,x)

velocity = 1;
dt = 1;
process_var = 1;

%Process Model
process_model = normpdf(t,9.7, process_var)
hold on
plot(t,process_model)

product_model = normpdf(t,98.94,process_var)
figure;plot(t,product_model,'--')
    
