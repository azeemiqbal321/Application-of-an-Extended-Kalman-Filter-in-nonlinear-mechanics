function update_out=updatem(x,measurement)
temp=x(2)+measurement(2);
update_out(1)=(x(1)*measurement(2)+x(2)*measurement(1))/temp;
update_out(2)=x(2)*measurement(2)/temp;
end