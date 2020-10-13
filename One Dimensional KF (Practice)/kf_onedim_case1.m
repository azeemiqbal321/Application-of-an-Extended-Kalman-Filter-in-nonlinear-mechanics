clear all
clc

t=[1:20];
process_var=1;
measurement_var=2;
velocity=1;
dt=1;
x=[1 401];  % Estimate of the position
dx=0;
pos = x(1);
var = x(2);
process=[velocity*dt process_var];
data = zeros(2,2);

for k=1:20
    z = x(1)+sqrt(process_var)*randn % Dog Simulation
    
    %predict
    dx = velocity*dt
    pos = pos + dx
    var = var + process_var

    %update
    pos  = (var*z + measurement_var*pos) / (var + measurement_var)
    var = (var * measurement_var) / (var + measurement_var)

    %save data
    data(k,1)=z;
    data(k,2)=pos;

end

plot(t,data(k,2),'o')
