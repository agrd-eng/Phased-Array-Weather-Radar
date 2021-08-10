function [X_hat,W] = CP(Y,steer_theta,Delta)
% CP Computed adaptive beamforming with gain control called Capon
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
%   [X_hat,W] = CP(Y,steer_theta,Delta)
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
%                 for a Uniform Linear Array (ULA) and lambda - system wavelength.[N,~] = size(Y);
[N,L] = size(Y);
M = numel(steer_theta); 

Ry = 1/L*(Y*Y'); % samples covariance matrix

for i = 1:M   
    S(:,i) = gen_a(N,Delta,steer_theta(i));
    W(:,i) = (Ry^(-1)*S(:,i))/((S(:,i)'*Ry^(-1)*S(:,i)));
end  

X_hat = W'*Y; 