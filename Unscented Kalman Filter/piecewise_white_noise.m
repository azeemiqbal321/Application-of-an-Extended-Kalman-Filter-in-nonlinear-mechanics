function output=piecewise_white_noise(ndim,variance,dt)

% We have defined F for the standard kinematic equations.
% USage: piecewise_white_noise(ndim=2 or 3,variance of noise,dt=time step)

if ndim==3
G=@(t)[t^2/2;t;1];
Gc=@(t)[t^2/2 t 1];
output=variance*G(dt)*Gc(dt);
else 

G=@(t)[t^2/2;t];
Gc=@(t)[t^2/2 t];
output=variance*G(dt)*Gc(dt);
end
end