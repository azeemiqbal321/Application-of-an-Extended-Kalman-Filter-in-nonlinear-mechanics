clc; clear all;
A = [0 1 0;-5 -2 5;0 0 1];
B = [0;3;0];
C = [1 0 0];
D = 0;
Ts = 0.05;
sys = ss(A,B,C,D,Ts,'StateName',{'Position' 'Velocity' 'Acceleration'},'InputName','Force');
Ob = obsv(A,C);
% Number of unobservable states
unob = length(A)-rank(Ob)
