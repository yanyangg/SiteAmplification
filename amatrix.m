function a=amatrix(alpha,beta,rho,p,w,z, idx)
% this is used to calculate a matrix(4*4)
% Reference: QUANTITATIVE SEISMOLOGY Page 166 (5.65)
% a=E*Lambda*E^-1
% Input:
% alpha - p velocity

% beta - s velocity
% rho - density
% p - ray parameter
% w - frequency
% z - z-z_ref textbook (5.56)
% idx - index, 0 for SH, 1 for P-SV
% Output:
% a - a matrix(4*4)

xi=sqrt(1/alpha^2-p^2);
eta=sqrt(1/beta^2-p^2);
miu = rho*beta^2;
if idx==1
    E=zeros(4,4);
    Lambda=zeros(4,4);
    % E matrix and Lambda matrix according to textbook (5.65)
    E=Ematrix(alpha,beta,rho,p,w,idx);

    Lambda(1,1)=exp(1i*w*xi*z);
    Lambda(2,2)=exp(1i*w*eta*z);
    Lambda(3,3)=exp(-1i*w*xi*z);
    Lambda(4,4)=exp(-1i*w*eta*z);

    a=E*Lambda*pinv(E);
else
    E=zeros(2,2);
    Lambda=zeros(2,2);
    % E matrix and Lambda matrix according to textbook (5.65)
    E=Ematrix(alpha,beta,rho,p,w,idx);

    Lambda(1,1)=exp(1i*w*eta*z);
    Lambda(2,2)=exp(-1i*w*eta*z);

    a=E*Lambda*pinv(E);
end
end