function xdot=duffing(t,x)

global gamma omega beta AMP rho

xdot(1)=x(2);
xdot(2)=-gamma*x(2)+omega^2*x(1)-beta*x(1)^3+AMP*cos(rho*t);

xdot=xdot';

end
