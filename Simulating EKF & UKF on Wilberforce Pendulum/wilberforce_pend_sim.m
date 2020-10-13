close all; clear all; clc;

% Assign model parameters
global omega epsilon mass inertia theta z
omega = 2.314;         % rad.s-1
epsilon = 9.27e-3;     % N
mass = 0.4905;         % kg
inertia = 1.39e-4;     % kg.m2
theta=2*pi;            % Initial angle
z=0;                   % Starting position

time = linspace(0,40,50000);
[t,y]=ode45(@wilberforce_func,time,[z 0 theta 0]);

figure(1)
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,3),'plot')
set(get(AX(1),'Ylabel'),'String','z (m)')
tix=get(AX(1),'ytick')';
set(gca,'yticklabel',num2str(tix,'%.2f'))
set(get(AX(2),'Ylabel'),'String','\theta (rad)')
legend 'Vertical Position (z)' 'Angular Position (\theta)'
xlabel('Time (s)')

O=[1 0 0 0;0 1 0 0;-omega^2 0 1/2/mass*epsilon 0;0 -omega^2 0 1/2/mass*epsilon]
rank(O)

% figure(2); plot(y(:,1),y(:,3))
% xlabel('z (m)'); ylabel('\theta (rad)')
% title('Phase Plot')
