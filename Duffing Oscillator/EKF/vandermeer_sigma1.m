function [chi,scalefactor,wm,wc]=vandermeer_sigma1(mu,sigma,alpha,beta,kappa)

%implementing van der Meer scaled sigma point calculation
%usage: [sigma,scalefector,wm,wc]=vandermeer_sigma1(mu,sigma,0.4,0.9,0.3)
%mu is column vector
%sigma is a square matrix

%generating the sigma points
ndim=length(mu); %dimension of state space
lambda=alpha^2*(ndim+kappa)-ndim;
chi=zeros(ndim,2*ndim+1); %pre-allocating space for the sigma points
chi(:,1)=mu; %Sigma points
%scalefactor=sqrtm((ndim+lambda)*sigma); %scale factor
scalefactor=chol((ndim+lambda)*sigma,'lower');

for k=2:ndim+1
    chi(:,k)=mu+scalefactor(:,k-1);
end

for k=ndim+2:2*ndim+1
    chi(:,k)=mu-scalefactor(:,k-ndim-1);
end

%generating the weights in means wm, and weights in covariances wc

wm=zeros(1,2*ndim+1);
wc=zeros(1,2*ndim+1);
wm(1)=lambda/(ndim+lambda);
wc(1)=wm(1)+1-alpha^2+beta;

for k=2:2*ndim+1
wm(k)=1/(2*(ndim+lambda));
wc(k)=wm(k);
end



end