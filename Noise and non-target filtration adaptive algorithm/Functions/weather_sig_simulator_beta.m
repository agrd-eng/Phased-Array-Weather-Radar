function [X_PS_nf, X_PS_n,M0_truth,M1_truth,M2_truth] = weather_sig_simulator_beta(M0,M1,M2,L,vmax,varn)
% WEATHER_SIG_SIMULATOR_BETA Simulates weather like time domain signals
%
% Description:
%       When considered that Doppler weather spectra, or any other point
%       target spectra can be parametrized as a Gaussian shape distribution
%       then it can be simulated. For this, the algorithm of Zrnic from [1]
%       is employed and implemented in this function. A gaussian
%       distribution can be characterized by three moments: M0, M1 and M2,
%       total power, mean and first standard deviation. When estimated
%       based on a Doppler spectra they are called: total power
%       reflectivity, mean Doppler velocity and spectral width.
%
% Usage:
%       [X_PS_nf, X_PS_n,M0_truth,M1_truth,M2_truth] = weather_sig_simulator_beta(M0,M1,M2,L,vmax,varn)
%
% Input:
%       M0       - vector, total power reflectivity;
%       M1       - vector, mean Doppler velocity;
%       M2       - vector, spectral width;
%       L        - scalar, number of Doppler bins;
%       vmax     - maximum unambiguous Doppler velocity;
%       varn     - noise variance when white noise.
% Output:
%       X_PS_nf  - vector, noise free time series signal;
%       X_PS_n   - vector, noisy time series signal;
%       M0_truth - vector, noisy total power reflectivity;
%       M1_truth - vector, noisy mean Doppler velocity;
%       M2_truth - vector, noisy spectral width.
%
% ========================================================================
% 18.08.2021 Helper added
% ========================================================================
M = numel(M0); % number of sources

for m = 1:M   
   
    ph = -pi+2*pi*rand(1,L);% uniformly distributed phase [-pi,pi]
    UD = rand(1,L);% uniformly distributed random variable [0,1]
    v = gauss_gen(vmax,L,M0(m),M1(m),M2(m));% v is in power
    
    X_PS_nf(m,:) = ((-v.*log(UD))).*exp(1i*ph);% noise free Power Spectrum
      
    X_PS_n(m,:) = ((-(v+varn(m)).*log(UD))).*exp(1i*ph);% Power spectrum

    [M0_truth(m,1),M1_truth(m,1),M2_truth(m,1)] = gauss_calc(abs(X_PS_n(m,:)),vmax,L);% Moments computation with noise

end

% [1] D. S. Zrnic,Simulation of Weatherlike Doppler Spectra and Signals,J.Appl.Meteorol.14, no.4, 619 (June 1975)
