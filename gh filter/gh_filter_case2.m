clear all
close all
clc

t=[1:12]; %No. of iterations
weights=randn(12,1); %Input Values
x_est = 5;
g = 6/10;
h = 6/10;
dx = 1;
dt = 1
data = zeros(4,4)

for k=1:12
   %prediction step
   x_pred = x_est + (dx*dt);
   dx=dx
   
   %update step
   residual = weights(k)-x_pred;
   dx = dx + h*(residual/dt)
   x_est = x_pred + g*residual; 
    
   %Saving Data
   %data(k,1)=t;
   data(k,2)=x_pred;
   data(k,3)=x_est;
   
end

%Plotting Data
plot(t,weights,'o')
hold on
plot(t,data(:,3),'-b','LineWidth',2)
legend('Measurements','Filter')
legend('Location','southeast')
title('The g-h filter implementation (case 2)')
xlabel('Time (s)')
ylabel('Weight (lbs)')

