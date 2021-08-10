
function f = gauss_gen(vmax,Ndoppler,M0,M1,M2)
%GAUSS_GEN Generates gaussian spectra based on input moments
%
%   Description:
%   Gaussian distribution is characterized by three moments, number of bins
%   and edges limits. This is used to generate a weatherlike Doppler
%   spectra which can be approximated with a Gaussian pdf shape.
%
%   Usage:
%   f = gauss_gen(vmax,Ndoppler,M0,M1,M2)
%   
%   Output:
%   f        - vector: 1 x Ndoppler, Gaussian probability distribution 
%              function values;
%
%   Input:
%
%   vmax     - scalar, maximum unambiguos Doppler velocity, |v_unamb|;
%   Ndoppler - scalar, number of Doppler bins;
%   M0       - scalar, input total power;
%   M1       - scalar, input mean Doppler velocity;
%   M2       - scalar, input spectral width, 1st std.
%
dv = 2*vmax/Ndoppler;% velocity resolution
x = -vmax:dv:vmax-dv;% velocity bins (contains zero frequency/velocity)

f = M0/(sqrt(2*pi)*M2)*exp(-(x-M1).^2/(2*M2^2));

end