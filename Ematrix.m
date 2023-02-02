function E=Ematrix(alpha,beta,rho,p,w, idx)
% this is used to calculate E matrix(4*4)
% Reference: QUANTITATIVE SEISMOLOGY Page 166 (5.65)
% Input:
% alpha - p velocity
% beta - s velocity
% rho - density
% p - ray parameter
% w - frequency
% idx - index, 0 for SH, 1 for P-SV
% Output:
% E - E matrix(4*4)
xi2=sqrt(1/alpha^2-p^2);
eta2=sqrt(1/beta^2-p^2);
miu2 = rho*beta^2;
if idx==1
    E=zeros(4,4);
    % E matrix according to textbook (5.65)
    E(1,1)=alpha*p;
    E(1,2)=beta*eta2;
    E(1,3)=alpha*p;
    E(1,4)=beta*eta2;

    E(2,1)=alpha*xi2;
    E(2,2)=-beta*p;
    E(2,3)=-alpha*xi2;
    E(2,4)=beta*p;

    E(3,1)=2*1i*w*rho*alpha*beta^2*p*xi2;
    E(3,2)=1i*w*rho*beta*(1-2*beta^2*p^2);
    E(3,3)=-2*1i*w*rho*alpha*beta^2*p*xi2;
    E(3,4)=-1i*w*rho*beta*(1-2*beta^2*p^2);

    E(4,1)=1i*w*rho*alpha*(1-2*beta^2*p^2);
    E(4,2)=-2*1i*w*rho*beta^3*p*eta2;
    E(4,3)=1i*w*rho*alpha*(1-2*beta^2*p^2);
    E(4,4)=-2*1i*w*rho*beta^3*p*eta2;
else
    E=zeros(2,2);
    E(1,1)=1;
    E(1,2)=1;
    E(2,1)=1i*w*miu2*eta2;
    E(2,2)=-1i*w*miu2*eta2;
end
end