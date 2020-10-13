data=load('data.lvm');

time=data(:,1)./60;
temp=data(:,2);

plot(time,temp)
xlabel('Time (s)');
ylabel('Temperature (°C)');
grid on;