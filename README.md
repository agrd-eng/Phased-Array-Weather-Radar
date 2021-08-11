# Main Repository

This repository contains personal developed algorithms for signal processing for radar systems. Different topics
are listed in different folders. In this README, each chapter contains a brief explanation for each folder.
More detailed explanations are found in each folder.

## Weatherlike and Phased Array Radar

This is a radar data simulator which first simulates Doppler spectra, of Gaussian shape, characterized by three moments:

1. Zero moment: total power reflectivity;
2. First moment: mean Doppler velocity [m/s];
3. Second moment: spectral width [m/s].

This is used to simulate time signals reflected by a weather object, but can be used for point targets as well (when spectral
width much smaller). These time signals are used as input for the data model of a phased array antenna, and the performances of
three different beamforming algorithms are tested:

1. Fourier
2. Capon
3. Recursive MMSE (Minimum Mean Square Error) with gain control [1]

Data model for the phased array antenna system is:
* Y = S*X+N
Where ground truth time signals are stored in matrix X, S is a matrix of array response vectors, a(theta), and N is i.i.d. Gaussian noise.
The role of the beamforming algorithms is to estimate X from Y, when Y is noisy.

## References
[1] E. Yoshikawa et al., “MMSE Beam Forming on Fast-Scanning Phased Array Weather Radar,” IEEE Trans. Geosci. Remote Sens., vol. 51, no. 2, pp. 3077–3088, 2013.