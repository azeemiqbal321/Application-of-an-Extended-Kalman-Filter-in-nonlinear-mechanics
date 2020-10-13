clear all; clc; 

syms x1 x2 x3 x4 x5 x6 x7 m I
x=[x1 x2 x3 x4 x5 x6 x7];
fx=[x2;-x5*x1-1/2/m*x7*x3;x4;-x6*x3-1/2/I*x7*x1;x5;x6;x7];
%y0=inline('x1');
Lie_0=x1;
% 
y1=inline('Lie_0');
Lie_1=[diff(y1(Lie_0),x1) diff(y1(Lie_0),x2) diff(y1(Lie_0),x3) diff(y1(Lie_0),x4) diff(y1(Lie_0),x5) diff(y1(Lie_0),x6) diff(y1(Lie_0),x7)]*fx;
% 
y2=inline('Lie_1');
Lie_2=[diff(y2(Lie_1),x1) diff(y2(Lie_1),x2) diff(y2(Lie_1),x3) diff(y2(Lie_1),x4) diff(y2(Lie_1),x5) diff(y2(Lie_1),x6) diff(y2(Lie_1),x7)]*fx;
% 
y3=inline('Lie_2');
Lie_3=[diff(y3(Lie_2),x1) diff(y3(Lie_2),x2) diff(y3(Lie_2),x3) diff(y3(Lie_2),x4) diff(y3(Lie_2),x5) diff(y3(Lie_2),x6) diff(y3(Lie_2),x7)]*fx;
% 
y4=inline('Lie_3');
Lie_4=[diff(y4(Lie_3),x1) diff(y4(Lie_3),x2) diff(y4(Lie_3),x3) diff(y4(Lie_3),x4) diff(y4(Lie_3),x5) diff(y4(Lie_3),x6) diff(y4(Lie_3),x7)]*fx;

y5=inline('Lie_4');
Lie_5=[diff(y5(Lie_4),x1) diff(y5(Lie_4),x2) diff(y5(Lie_4),x3) diff(y5(Lie_4),x4) diff(y5(Lie_4),x5) diff(y5(Lie_4),x6) diff(y5(Lie_4),x7)]*fx;

y6=inline('Lie_5');
Lie_6=[diff(y6(Lie_5),x1) diff(y6(Lie_5),x2) diff(y6(Lie_5),x3) diff(y6(Lie_5),x4) diff(y6(Lie_5),x5) diff(y6(Lie_5),x6) diff(y6(Lie_5),x7)]*fx;

y7=inline('Lie_6');
Lie_7=[diff(y7(Lie_6),x1) diff(y7(Lie_6),x2) diff(y7(Lie_6),x3) diff(y7(Lie_6),x4) diff(y7(Lie_6),x5) diff(y7(Lie_6),x6) diff(y7(Lie_6),x7)]*fx;

phi=[Lie_0;Lie_1;Lie_2;Lie_3;Lie_4;Lie_5;Lie_6;Lie_7];

O=[diff(Lie_0,x1) diff(Lie_0,x2) diff(Lie_0,x3) diff(Lie_0,x4) diff(Lie_0,x5) diff(Lie_0,x6) diff(Lie_0,x7);
   diff(Lie_1,x1) diff(Lie_1,x2) diff(Lie_1,x3) diff(Lie_1,x4) diff(Lie_1,x5) diff(Lie_1,x6) diff(Lie_1,x7);
   diff(Lie_2,x1) diff(Lie_2,x2) diff(Lie_2,x3) diff(Lie_2,x4) diff(Lie_2,x5) diff(Lie_2,x6) diff(Lie_2,x7);
   diff(Lie_3,x1) diff(Lie_3,x2) diff(Lie_3,x3) diff(Lie_3,x4) diff(Lie_3,x5) diff(Lie_3,x6) diff(Lie_3,x7);
   diff(Lie_4,x1) diff(Lie_4,x2) diff(Lie_4,x3) diff(Lie_4,x4) diff(Lie_4,x5) diff(Lie_4,x6) diff(Lie_4,x7);
   diff(Lie_5,x1) diff(Lie_5,x2) diff(Lie_5,x3) diff(Lie_5,x4) diff(Lie_5,x5) diff(Lie_5,x6) diff(Lie_5,x7);
   diff(Lie_6,x1) diff(Lie_6,x2) diff(Lie_6,x3) diff(Lie_6,x4) diff(Lie_6,x5) diff(Lie_6,x6) diff(Lie_6,x7);
   diff(Lie_7,x1) diff(Lie_7,x2) diff(Lie_7,x3) diff(Lie_7,x4) diff(Lie_7,x5) diff(Lie_7,x6) diff(Lie_7,x7)]



