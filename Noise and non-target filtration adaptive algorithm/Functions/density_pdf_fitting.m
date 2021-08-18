function rmse = density_pdf_fitting(density_in)
%DENSITY_PDF_FITTING Compute rmse between estimated and U-quadratic pdf
%
% Description
%                     Computes RMSE value between estimated pdf
%                     given by density vector returned by binary_scan_win.m 
%                     function: density_in,, and pdf for U-quadratic
%                     distribution.
%                     Parametric U-quadratic pdf is computed considering
%                     that density range is between [0,1]. In order to
%                     computed non-parametric pdf and to estimate rmse
%                     value, which is a scalar, we use a kernel pdf estimation 
%                     technique with adaptive width.
%                     This can be found under MATLAB functions as fitdist.m
%                     and then pdf.m to generate the pdf. These two
%                     functions are packed under function called
%                     Kernel_pdf.m which can be found in MAX3D Toolbox. 
% Usage:
%       rmse = density_pdf_fitting(density_in)
% Input:
%       density_in  - vector, pixels density values computed for each position of a
%                     fixed-size sliding window, in a raw binary image.
% Output
%       rmse        - root mean square error - scalar value computed
%                     between estimated pdf using density_in using a Kernel
%                     estimator (Gaussian kernel), and U-quadratic pdf.
%==========================================================================
% v.1.0 - AG, 2021
% 02.06.2021, OK - TODO - add more descriptions
% August - AG, description added
%==========================================================================
a = 0;
b = 1;

x = linspace(a,b,100);

beta = (a+b)/2;
alpha = 12/(b-a)^3;
u_quadr = alpha*(x-beta).^2;

y_out = Kernel_pdf(density_in,x);

rmse = sqrt(mean((u_quadr-y_out).^2));

end

