function xdot=wilberforce_func2(t,x)

global omega gamma mass inertia epsilon  

xdot(1)=x(2);
xdot(2)=(-omega^2*x(1)-0.5/mass*epsilon*x(3))*t;
xdot(3)=x(4);
xdot(4)=(-gamma^2*x(3)-0.5/inertia*epsilon*x(1))*t;

xdot=xdot';

end
