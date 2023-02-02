function u = site_resp_body(thk,dns,vp,vs,freq,rayp,flag)
% This function is used to compute site amplification term for body wave
% using propagation matrix (Aki & Richards, 2002).
% Usage: propagation_mat, amatrix.m, Ematrix.m
% Input:
% thk - layer thickness, N-1 * 1
% dns - density, N * 1
% vp - P velocity, N * 1
% vs - S velocity, N * 1
% freq - a vector of desired frequency, nf * 1
% rayp - ray parameter for incident wave.
% flag - 'P' for P wave, 'SV' for SV wave, 'SH' for SH wave.
% Output:
% u - ground displacement, if flag='R', u is nf * 2, where 1st column is
% the radial component, 2nd column is vertical component. If flag='L', u 
% is nf * 1, the transverse component.
nf = length(freq);
w = freq*2*pi;
if strcmp(flag,'P')
    ur = zeros(nf,1);
    uz = zeros(nf,1);
    for i=1:nf
        M=propagation_mat(vp,vs,dns,rayp,w(i),thk, 1); 
        A=M(1:2,1:2);
        B=M(1:2,3:4);
        C=M(3:4,1:2);
        D=M(3:4,3:4);
        R=A-B/D*C;
        ur(i)=R(1,1);
        uz(i)=R(2,1);
    end
    u = [ur, uz];
elseif strcmp(flag,'SV')
    ur = zeros(nf,1);
    uz = zeros(nf,1);
    for i=1:nf
        M=propagation_mat(vp,vs,dns,rayp,w(i),thk, 1); 
        A=M(1:2,1:2);
        B=M(1:2,3:4);
        C=M(3:4,1:2);
        D=M(3:4,3:4);
        R=A-B*pinv(D)*C;
        ur(i)=R(1,2);
        uz(i)=R(2,2);
    end
    u = [ur, uz];
elseif strcmp(flag,'SH')
    ut = zeros(nf,1);
    for i=1:nf
        M=propagation_mat(vp,vs,dns,rayp,w(i),thk, 0); 
        ut(i)=M(1,1)-M(1,2)*M(2,1)/M(2,2);
    end
    u = ut;
else
    error('flag should be P, SV or SH!')
end