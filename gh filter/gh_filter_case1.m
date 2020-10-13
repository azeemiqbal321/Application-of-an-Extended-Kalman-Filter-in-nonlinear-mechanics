clear all
close all
clc

t=[1:12]; %No. of iterations
weights=[158.0 164.2 160.3 159.9 162.1 164.6 169.6 167.4 166.4 171.0 171.2 172.6]; %Input Values
x_est = 160;
g = 6/10;
h = 2/3;
dx = 1;
dt = 1
data = zeros(4,4)

for k=1:12
   %prediction step
   x_pred = x_est + (dx*dt);
   
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
title('The g-h filter implementation (case 1)')
xlabel('Time (s)')
ylabel('Weight (lbs)')

