function A = gen_a(M,Delta,theta)
%GEN_A Generates array response vector, a(theta(i))
%
%   Description
%   Model used for phased array antenna is:
%   X = A*S+N
%   X: M x N - received data, M number of antenna elements and N - number
%              of samples
%   A: M x k - arrays response vectors on column
%   S: k x N - signals reflected by each target at k^{th} time sample, row
%   N: M x N - Gaussian additive noise
%
%  Usage:
%  A = gen_a(M,Delta,theta)
%
%  Output:
%  A      - matrix, M x N, where M is the number of antenna elements and N is
%           the number of time samples. Each column represent a(theta(k)),
%           for each k^{th} time sample.
%
%  Input:
%  M      - scalar, number of antenna elements.
%  Delta  - scalar, ration between d - distance between antenna elements
%           for a Uniform Linear Array (ULA) and lambda - system wavelength. 
%  theta  - steering vector, 1 x k, where k is the number of steering
%           angles.

A = exp((0:M-1).'*1i*2*pi*Delta*sind(theta));
end

