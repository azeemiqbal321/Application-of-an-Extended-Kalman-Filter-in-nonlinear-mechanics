clc
close all
clear all

%% Simulating Duffing Oscillator
global gamma omega beta F rho

gamma=0.3;  %Damping Coefficient
omega=-1;    %Stiffness Coefficient
beta=1;  %Stiffness Coefficient
F=0.5;  %Amplitude of driving force
rho=1.25;  %Frequency of driving force

[t x]=ode45(@duffing,0:2*pi/rho/100:4000,[1 0]);

dst=t(1:6000);  %time
yn=x(1:6000,1); %position
yv=x(1:6000,2); %velocity

%For Poincare Section
for i=200:100:78000
     n=(i-100)/100;
     x1(n)=x(i,1);
     x2(n)=x(i,2);
end

%% Plotting
subplot(2,2,[1 2])
plot(dst,yn,'-')
axis([0 300,-2 2]);
tix=get(gca,'ytick')';
set(gca,'fontsize',10)
set(gca,'yticklabel',num2str(tix,'%.1f'))
xlabel('Time (s)')
ylabel('$\theta (rad)$','interpreter','latex')
legend('True Position','Position with noise','Linear Kalman Filter')

%Poincare Section
subplot(2,2,3)
plot(x1(:),x2(:),'*')
axis([-1.5 1.3,-0.5 1.0]);
tix=get(gca,'ytick')';
set(gca,'fontsize',10)
set(gca,'yticklabel',num2str(tix,'%.1f'))
xlabel('\theta (rad)')
ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex')

%Phase space
subplot(2,2,4)
plot(yn,yv,'b')
axis([-1.8 1.8,-1.2 1.2]);
tix=get(gca,'ytick')';
set(gca,'fontsize',10)
set(gca,'yticklabel',num2str(tix,'%.1f'))
xlabel('\theta (rad)')
ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex')