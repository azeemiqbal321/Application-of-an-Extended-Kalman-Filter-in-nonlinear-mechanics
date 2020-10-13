clear all; clc; close all;

omega = 2.314;         % rad.s-1
epsilon = 9.27e-3;     % N
mass = 0.4905;         % kg
inertia = 1.39e-4;     % kg.m2
theta=0;            % Initial angle
z=0;                   % Starting position

A=[1 0 0 0;
   0 1 0 0;
   -omega^2 0 0.5/2*mass*epsilon 0 ;
   0 -omega^2 0 0.5/2*mass*epsilon]

rank(A)


