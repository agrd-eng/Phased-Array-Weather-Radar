function [X_hat,iter,W] = cMMSE(Y,steer_theta,Delta,varn,th)
% cMMSE Computed adaptive re-iterative beamforming method
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
%   [X_hat,iter,W] = cMMSE(Y,steer_theta,Delta,varn,th)
%
%   Output:
%   X_hat       - signals matrix, estimated reflected signals
%   iter        - scalar, returns the number of iterations necessary to
%                 exceed the threshold
%   W           - beamforming matrix
%
%   Input:
%   Y           - data matrix, samples received by the phased array antenna system
%   steer_theta - steering vector, 1 x M, M here is the total number of
%                 steering angles
%   Delta       - scalar, ration between d - distance between antenna elements
%                 for a Uniform Linear Array (ULA) and lambda - system wavelength.[N,~] = size(Y);
%   varn        - scalar, noise variance, used to compute noise covariance
%                 matrix, when noise is i.i.d.
%   th          - scalar, fixed value, 0.001 proposed in [1]

[N,L] = size(Y);
M = numel(steer_theta); 
Rv = varn*eye(N); % noise correlation matrix

for i = 1:M   
    S(:,i) = gen_a(N,Delta,steer_theta(i));
end  
W = 1/N*S; % MxN :FR beamformer

iter = 0;
delta_th = 1;
Rx = zeros(M,M);
temp = zeros(M,M);

while delta_th > th
    if iter > 1
        X_hat_temp = X_hat;
    end
    X_hat = W'*Y; % MxN*NxL ->X_hat dim: MxL    
    Rx = 1/L.*X_hat*X_hat'.*eye(M);    
    R = S*Rx*S'+Rv;   
    clear W;
    for ii = 1:M
        W(:,ii) = (R^(-1)*S(:,ii))/(S(:,ii)'*R^(-1)*S(:,ii));
    end    
    if iter > 1
        for iii = 1:M
            temp2(iii) = norm(X_hat(iii,:)-X_hat_temp(iii,:),2)^2/norm(X_hat_temp(iii,:))^2;
        end
        delta_th = 1/M*sum(temp2);
    end
    
    iter = iter+1; 
end

% [1] E. Yoshikawa et al., “MMSE Beam Forming on Fast-Scanning 
% Phased Array Weather Radar,” IEEE Trans. Geosci. Remote Sens., 
% vol. 51, no. 2, pp. 3077–3088, 2013.
