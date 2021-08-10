function [X_TD,X_PS_n,M0_truth,M1_truth,M2_truth] = weatherlike_spectra_signals(M0,M1,M2,L,vmax,varn)
%WEATHERLIKE_SPECTRA_SIGNALS Simulate weatherlike Doppler spectra and signals
%
%   Description:
%   To study the performances of different algorithms used for radar weather
%   applications, e.g. beamforming, in case of phased array radar, weather targets
%   detection and noise or/and clutter cancelation algorithms etc., a
%   simulation technique is required. Within this simulator Doppler spectra
%   is generated as a Gaussian shape, which requires three moments: total
%   power reflectivity, mean Doppler velocity and spectral width. Moreover,
%   additive noise power can be introduced in the spectra. Formula and
%   algorithm used is taken from [1].
% 
%   Usage:
%   [X_TD,X_PS_n,M0_truth,M1_truth,M2_truth] = weatherlike_spectra_signals(M0,M1,M2,L,vmax,varn)
%
%   Output:
%   X_TD     - matrix: K x L, matrix containing weather like signals - time
%              domain, for L Doppler bins and K is either the number of range bin 
%              or elevation bins, whatever user wants to generate
%   X_PS     - matrix: K x L, matrix containing weather like Doppler spectra
%              - frequency domain;
%   M0_truth - matrix: K x 2, that contains total power (zero moment)
%              values for weather like Doppler spectra with additive noise 
%              for the 1st row and noise free for the 2nd one;
%   M1_truth - matrix: K x 2, that contains mean Doppler velocity (1st moment)
%              values for weather like Doppler spectra with additive noise 
%              for the 1st row and noise free for the 2nd one;
%   M2_truth - matrix: K x 2, that contains spectral width (2nd moment)
%              values for weather like Doppler spectra with additive noise 
%              for the 1st row and noise free for the 2nd one.
%  
%   Input:
%   M0       - vector, 1 x K: total power reflectivity;
%   M1       - vector, 1 x K: mean Doppler velocity;
%   M2       - vector, 1 x K: spectral width (1st standard deviation)
%   L        - number of Doppler bins
%   vmax     - maximum unambiguos Doppler velocity, |v_unamb|
%   varn     - noise variance

K = numel(M0); % number of range bins/elevation bins

for m = 1:K   
   
    ph = -pi+2*pi*rand(1,L);% uniformly distributed phase [-pi,pi]
    UD = rand(1,L);% uniformly distributed random variable [0,1]
    v = gauss_gen(vmax,L,M0(m),M1(m),M2(m));% v is returned for power values
    
    X_PS_nf(m,:) = ((-(v).*log(UD))).*exp(1i*ph);% power spectrum [1]
    X_TD(m,:) = ifft(fftshift((L.*X_PS_nf(m,:)).^(1/2))); % power spectrum -> amplt spectrum(denormalized) -> signal time domain
   
    X_PS_n(m,:) = ((-(v+varn).*log(UD))).*exp(1i*ph);% power spectrum [1]

   
    [M0_truth(m,1),M1_truth(m,1),M2_truth(m,1)] = gauss_calc(abs(X_PS_n(m,:)),vmax);% Moments computation with noise
    [M0_truth(m,2),M1_truth(m,2),M2_truth(m,2)] = gauss_calc(abs(X_PS_nf(m,:)),vmax);% noise-free Moments computation

end
% [1] D. S. Zrnic,Simulation of Weatherlike Doppler Spectra and Signals,J.Appl.Meteorol.14, no.4, 619 (June 1975)

end

