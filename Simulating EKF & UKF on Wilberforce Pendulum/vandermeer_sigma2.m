function [chi,scalefactor,wm,wc]=vandermeer_sigma2(muu,sigma,alpha,beta,kappa)

%implementing van der Meer scaled sigma point calculation
%usage: [sigma,scalefector,wm,wc]=vandermeer_sigma1(mu,sigma,0.4,0.9,0.3)
%mu is column vector
%sigma is a square matrix

%generating the sigma points
ndim2=length(muu); %dimension of state space
lambda=alpha^2*(ndim2+kappa)-ndim2;
chi=zeros(ndim2,2*ndim2+1); %pre-allocating space for the sigma points
chi(:,1)=muu;

scalefactor=real(sqrtm((ndim2+lambda)*sigma));
%scalefactor=chol((ndim2+lambda)*sigma);

for k=2:ndim2+1
    chi(:,k)=muu+scalefactor(:,k-1);
end

for k=ndim2+2:2*ndim2+1
    chi(:,k)=muu-scalefactor(:,k-ndim2-1);
end

%generating the weights in means wm, and weights in covariances wc

wm=zeros(1,2*ndim2+1);
wc=zeros(1,2*ndim2+1);
wm(1)=lambda/(ndim2+lambda);
wc(1)=wm(1)+1-alpha^2+beta;

for k=2:2*ndim2+1
wm(k)=1/(2*(ndim2+lambda));
wc(k)=wm(k);
end

 %figure; error_ellipse(sigma,mu); hold on; scatter(chi(1,:),chi(2,:));

end