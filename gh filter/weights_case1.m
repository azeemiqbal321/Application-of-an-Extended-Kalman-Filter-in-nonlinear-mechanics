clear all
close all
clc

t=[1:12]; %No. of iterations
weights=[158.0 164.2 160.3 159.9 162.1 164.6 169.6 167.4 166.4 171.0 171.2 172.6]; %Input Values
Initial_guess=160;
weight=Initial_guess; %Initial guess
gain_rate=1;
time_step = 1; %day
scale_factor = 4/10;
data = ones(4,1);

for k=1:12
    
    prediction = weight + gain_rate*time_step; 
    weight = prediction + (scale_factor * (weights(k) - prediction));
    
    %Saving Data
    data(1,1)=Initial_guess;
    data(k+1:12,1)=weight;
    data(k,2)=prediction;
    data(k,3)=weight;
    
end

%Plotting Data
plot(t,weights,'o')
hold on
plot(t,data(:,2),'-b','LineWidth',2)
plot(t,data(:,3),'--r')
legend('Measurements','Estimations','Predictions')
legend('Location','southeast')
title('The weight case 1')
xlabel('Day')
ylabel('Weight (lbs)')

    