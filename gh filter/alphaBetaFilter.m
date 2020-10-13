function data = alphaBetaFilter(samples,m, x_est, dx,dt, g, h)
data=zeros(1,1);
for k=1:samples
   %prediction step
   x_pred = x_est + (dx*dt);  
   
   %update step
   residual = m(k)-x_pred;
   dx = dx + h*(residual/dt);
   x_est = x_pred + g*residual; 
   
   %save data
   data(k,1)=x_est;
end
end
