function [X_hat,W] = FR(Y,steer_theta,Delta)
% FR Computed classic beamforming, called Fourier or MF
%
%   Description:
%   Model:
%   Y = S*X+N 
%   Y dim: NxL
%   S dim: NxM
%   X dim: MxL
%   N: number of antenna elements
%   L: number of samples
%   M: number of steering angles
%
%   Usage:
%   [X_hat,W] = FR(Y,steer_theta,Delta)
%
%   Output:
%   X_hat       - signals matrix, estimated reflected signals
%   W           - beamforming matrix
%
%   Input:
%   Y           - data matrix, samples received by the phased array antenna system
%   steer_theta - steering vector, 1 x M, M here is the total number of
%                 steering angles
%   Delta       - scalar, ration between d - distance between antenna elements
%                 for a Uniform Linear Array (ULA) and lambda - system wavelength.
[N,~] = size(Y);
S = gen_a(N,Delta,steer_theta);

W = 1/N*S; % MxN :FR beamformer
X_hat = W'*Y; % MxN*NxL ->X_hat dim: MxL

end


