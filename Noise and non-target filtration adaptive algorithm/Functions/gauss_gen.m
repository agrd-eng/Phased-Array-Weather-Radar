function f = gauss_gen(vmax,Ndoppler,M0,M1,M2)
%GAUSS_GEN Generates Gaussian distribution shape
% Description:
%   To compute Doppler spectrum moments we assume weather Doppler spectra
%   follows a Gaussian distribution. Thus zero moment, M0 is total power
%   reflectivity, first moment, M1, is mean Doppler velocity, and second
%   moment, M2, is first standard deviation (aka spectral width). Using
%   this inpus parameters, number of bins (= NDoppler) and edges |v_{max}|,
%   one can generate any Gaussian distribution.
%Usage:
%        f = gauss_gen(vmax,Ndoppler,M0,M1,M2)
%Input
%   vmax     - scalar, maximum unambigous Doppler velocity;
%   Ndoppler - scalar, number of Doppler bins;
%   M0,M1 and M2 - scalar values, represents total power, mean Doppler
%   velocity and spectral width respectively.
%
%Output
%   f        - vector, resulted values for the Gaussian distribution shape.
%==========================================================================
% v.1.0 - AG, 2021
% August 2021 - Helper added
%==========================================================================
dv = 2*vmax/(Ndoppler);% velocity resolution
x = -vmax:dv:vmax-dv;% velocity bins

f = M0/(sqrt(2*pi)*M2)*exp(-(x-M1).^2/(2*M2^2));
