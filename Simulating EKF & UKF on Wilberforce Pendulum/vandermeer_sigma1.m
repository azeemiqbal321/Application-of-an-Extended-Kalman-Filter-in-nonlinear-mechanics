function [chi,scalefactor,wm,wc]=vandermeer_sigma1(mu,sigma,alpha,beta,kappa)

%implementing van der Meer scaled sigma point calculation
%usage: [sigma,scalefector,wm,wc]=vandermeer_sigma1(mu,sigma,0.4,0.9,0.3)
%mu is column vector
%sigma is a square matrix

%generating the sigma points
ndim=length(mu); %dimension of state space
lambda=alpha^2*(ndim+kappa)-ndim;
chi=zeros(ndim,2*ndim); %pre-allocating space for the sigma points
chi(:,1)=mu;
scalefactor=sqrtm((ndim+lambda)*sigma); %scale factor

for k=2:ndim
    chi(:,k)=mu+scalefactor(:,k);
end

for k=ndim+1:2*ndim
    chi(:,k)=mu-scalefactor(:,k-ndim);
end

%generating the weights in means wm, and weights in covariances wc

wm=zeros(1,2*ndim);
wc=zeros(1,2*ndim);
wm(1)=lambda/(ndim+lambda);
wc(1)=wm(1)+1-alpha^2+beta;

for k=2:2*ndim
wm(k)=1/(2*(ndim+lambda));
wc(k)=wm(k);
end

 %figure; error_ellipse(sigma,mu); hold on; scatter(chi(1,:),chi(2,:));

end