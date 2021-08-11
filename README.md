# Main Repository

This repository contains personal developed algorithms for signal processing for radar systems. Different topics
are listed in different folders. In this README, each chapter contains a brief explanation for each folder.
More detailed explanations can be found in each folder.

## Weatherlike and Phased Array Radar

This is a phased array radar data simulator for both extended (weather) and point targets, which can be tested
for fast rotating radar. The final goal is to test the performances of three different beamforming (BF) methods:

1. Fourier; [1]
2. Capon; [1]
3. Recursive MMSE (Minimum Mean Square Error) with gain control [1];

w.r.t. weather objects, when three main criteria for the BF methods are researched:

1. Robusted to a limited number of Doppler bins (i.e. fast rotating radar);
2. Accurate estimation of amplitude and phase information;
3. Performances for sidelobes suppression.

This goal is achieved by generating desired Doppler weather spectra, characterized by three statistical moments:
1. Zero moment: total power reflectivity;
2. First moment: mean Doppler velocity [m/s];
3. Second moment: spectral width [m/s];
where Doppler weather spectra is assumed to follow a Gaussian distribution shape. Next, Doppler data is transformed 
(using ifft) in time domain and used as input for the the data model used for a phased array antenna:
* Y = S*X+N
,where time domain signals (ground truth) are stored in matrix X, S is a matrix of array response vectors, a(theta(k)), 
and N is i.i.d. Gaussian noise. The role of the beamforming algorithms is to estimate X from Y, when Y is noisy.

NOTE: weather object represents an object spread over multiple elevation angles with wide spectra width. However, point
targets can be simulated as well, when elevation spread and spectral width are much narrower.

Finally, \hat{X} for each BF method is processed for Doppler domain and statistical moments are computed, such that,
the performances of the BF methods are compared between them and with the ground truth moments values.


## References
[1] E. Yoshikawa et al., “MMSE Beam Forming on Fast-Scanning Phased Array Weather Radar,” IEEE Trans. Geosci. Remote Sens., vol. 51, no. 2, pp. 3077–3088, 2013.