function output=continuous_white_noise(spectral_density,time_step)

% We have defined F for the wilberforce linear and angular movement
% [position,velocity,angle,angular velocity]^T as the state function


F=@(t)[1 t 0 0;0 1 0 0;0 0 1 t;0 0 0 1];
Fc=@(t)[1 0 0 0;t 1 0 0;0 0 1 0;0 0 t 1];
Qc=[0 0 0 0;0 1 0 0;0 0 0 0;0 0 0 1];
output=integral(@(t) spectral_density*F(t)*Qc*Fc(t),0,time_step,'ArrayValued',true);
end