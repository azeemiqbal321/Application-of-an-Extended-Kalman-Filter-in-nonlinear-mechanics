function xdot=duffing(t,x)

global gamma omega epsilon AMP OMEG

xdot=zeros(2,1);
xdot(1)=x(2);
xdot(2)=-gamma*x(2)+omega^2*x(1)-epsilon*x(1)^3+AMP*cos(OMEG*t);

end
