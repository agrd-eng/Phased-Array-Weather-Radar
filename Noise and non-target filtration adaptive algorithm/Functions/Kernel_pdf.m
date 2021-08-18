
function y_out = Kernel_pdf(input,x,BW)
%KERNEL_PDF Computes continuous non-parametric pdf using Kernel estimator technique. 
%
% Description:
%            Kernel_pdf(input,x,Bw) function pack together two MATLAB functions
%            PD = fitdist(input,BW) and pdf(PD,x), where the former fits 
%            the input values to a specific distribution, in this case a 
%            kernel and by default is 'gaussian', and  width: BW. If BW is not
%		     set then is adaptively computed. This function, "fitdist" returns PD,
%			 which plugged in latter function returns an array
%            of values of the pdf, for the pdf specified in PD evaluated at
%            the values in x.
% Usage:
% 			y_out = Kernel_pdf(input,x,BW)
% Input:
%			input    - vector, contains data samples;
%			x     	 - vector, contains data axis;
%			BW 	     - scalar, represent width, which is the number of samples
%                      over which each kernel distribution is computed.  
% Output:
%			y_out    - vector, final resulted distribution, using bassically 
%					   kernel estimation method, which by default for fitdist 
%					   uses Gaussian shape.
%==========================================================================
% v.1.0 - AG, 2021
% August 2021 - Helper Added
%==========================================================================
if nargin == 3
    pd_input = fitdist(input,'Kernel','BandWidth',BW);
else
    pd_input = fitdist(input,'Kernel');
end

y_out = pdf(pd_input,x);

end