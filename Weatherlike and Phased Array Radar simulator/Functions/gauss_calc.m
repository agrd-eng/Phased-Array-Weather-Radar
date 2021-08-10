function [M0t,M1t,M2t] = gauss_calc(v,vmax)
%GAUSS_CALC Computes Gaussian pdf moments based on the input pdf
%
%   Description:
%   Gaussian distribution is characterized by three moments:
%   M0: total power
%   M1: mean value
%   M2: 1st standard deviation
%   This is used to generate a weatherlike Doppler spectra for data 
%   received by a radar. Weather spectra can be approximated with a 
%   Gaussian pdf shape, we talk here about spectra moments:
%   M0: total power reflectivity
%   M1: mean Doppler velocity (considered redial wind velocity, v_r [m/s])
%   M2: spectral width [m/s]
%
%   Usage:
%   [M0t,M1t,M2t] = gauss_calc(v,vmax)
%   
%   Output:
%   M0t      - scalar, total power reflectivity, [linear scale]
%   M1t      - scalar, mean Doppler velocity, [m/s]
%   M2t      - scalar, spectral width, [m/s]
%
%   Input:
%   v        - vector: 1 x N, where N is the number of Doppler bins;
%   vmax     - scalar, maximum unambiguos Doppler velocity, |v_unamb|.
%

Ndoppler = numel(v);
dv = 2*vmax/Ndoppler;
x = -vmax:dv:vmax-dv;

M0t = sum(v*dv);
M1t = (1/M0t)*sum(x.*abs(v).*dv);
M2t = sqrt((1/(M0t)*sum((x-M1t).^2.*abs(v).*dv)));

end