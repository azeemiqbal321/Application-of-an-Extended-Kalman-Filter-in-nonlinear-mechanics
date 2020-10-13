function xdot=damped_ho(t,x)

global kspring bspring mass

xdot=zeros(2,1);

xdot(1)=x(2);
xdot(2)=-(kspring/mass)*x(1)-(bspring/mass)*x(2);
end
