function M=propagation_mat(alpha,beta,rho,p,w,thi, idx)
% this is used to calculate a matrix (2*2 for SH, 4*4 for P-SV) M
% where df/dz =  Mf
% Usage: amatrix.m, Ematrix.m
% Reference: QUANTITATIVE SEISMOLOGY Page 166 (5.65)
% Input:
% alpha - p velocity
% beta - s velocity
% rho - density
% p - ray parameter
% w - frequency
% thi - z-z_ref textbook (5.56)
% idx - index, 0 for SH, 1 for P-SV
% Output:
% M - a matrix (2*2 for SH, 4*4 for P-SV)
nlayer = length(alpha);
if idx==0
    a_prod=eye(2);
else
    a_prod=eye(4);
end
for j=1:nlayer-1
    a=amatrix(alpha(j),beta(j),rho(j),p,w,thi(j),idx);
    a_prod=a_prod*a;
end
E=Ematrix(alpha(end),beta(end),rho(end),p,w,idx);
M=a_prod*E;

