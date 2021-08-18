function [M0t,M1t,M2t] = gauss_calc(v,vmax,Ndoppler)
%GAUSS_CALC Calculates Doppler spectrum moments
% Description:
%   To compute Doppler spectrum moments we assume weather Doppler spectra
%   follows a Gaussian distribution. Thus zero moment, M0 is total power
%   reflectivity, first moment, M1, is mean Doppler velocity, and second
%   moment, M2, is first standard deviation (aka spectral width).
%Usage:
%        [M0t,M1t,M2t] = gauss_calc(v,vmax,Ndoppler)
%Input
%   v        - Doppler velocity grid;
%   vmax     - maximum unambiguous Doppler velocity;
%   Ndoppler - number of Doppler bins.
%Output
%   M0t - zero moment of Gausian Doppler spectrum (power)
%   M1t - first moment of Gausian Doppler spectrum (mean velocity)
%   M2t - second moment of Gausian Doppler spectrum (power)
%==========================================================================
% v.1.0 - AG, 2021
% 02.06.2021, OK - TODO - add more descriptions
%==========================================================================
dv = 2*vmax/Ndoppler;
x = -vmax:dv:vmax-dv;

M0t = sum(v*dv);
M1t = (1/M0t)*sum(x.*abs(v).*dv);
M2t = sqrt((1/(M0t)*sum((x-M1t).^2.*abs(v).*dv)));

end