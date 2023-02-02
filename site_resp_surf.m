function u = site_resp_surf(thk,dns,vp,vs,freq,flag)
% This function is used to compute site amplification term for surface 
% wave using Bodwen & Tsai (2016).
% Usage: the mat_disperse package to calculate the surface wave
% eigenfunctions.
% Input:
% thk - layer thickness, N-1 * 1
% dns - density, N * 1
% vp - P velocity, N * 1
% vs - S velocity, N * 1
% freq - a vector of desired frequency, M * 1
% flag - 'R' for Rayleigh wave, 'L' for Love wave.
% Output:
% u - ground displacement, if flag='R', u is M * 2, where 1st column is
% the radial component, 2nd column is vertical component. If flag='L', u 
% is M * 1, the transverse component
M = length(freq);
if strcmp(flag,'R')
    [vr,dvrvs,dvrrho,dvrvp,r,z,zdvrvs,zdvrrho,zdvrvp,I1,I2,I3,U] = mat_disperse(thk,dns,vp,vs,freq,'R');
    ur = abs(squeeze(r(:,1,1,1)));
    uz = abs(squeeze(r(:,1,1,2)));
    uzz = (abs(U).*abs(I1)).^(-0.5).*uz;
    urr = (abs(U).*abs(I1)).^(-0.5).*ur;
    u = [reshape(urr, M, 1), reshape(uzz, M, 1)];
elseif strcmp(flag,'L')
    [vr,dvrvs,dvrrho,r,z,zdvrvs,zdvrrho,I1,I2,I3,U] = mat_disperse(thk,dns,vp,vs,freq,'L');
    ut = abs(squeeze(r(:,1,1,1)));
    utt = (abs(U).*abs(I1)).^(-0.5).*ut;
    u = reshape(utt, M, 1);
else
    error('flag should be R or L!')
end