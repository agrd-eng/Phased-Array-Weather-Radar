function density_out = binary_scan_win(binary_data_in,win_s)
%BINARY_SCAN_WIN Computes pixel density in a binary image using a window
%
% Description:
% Computes pixel density within a fixed-size sliding window.  Pixel density 
% is ratio of pixels of one and window area.
%
% Usage: 
%         density_out = binary_scan_win(binary_data_in,win_s)
% Input:
%         binary_data_in - 2D binary image, [range] x [Doppler bins].
%         win_s          - window size, scalar.
% Output:
%         density_out    - pixel density values for each position of the
%                          sliding window, vector.
%==========================================================================
% v.1.0 - AG, 2021
% 22.07.2021, AG - Help info
%==========================================================================

sd = size(binary_data_in);

win_sv = win_s;
win_sh = win_s;

while mod(sd(1),win_sv)~=0
    win_sv = win_sv+1;
end

while mod(sd(2),win_sh)~=0
    win_sh = win_sh+1;
end

max_den = win_sh*win_sv;

sh = sd(2)/win_sh;
sv = sd(1)/win_sv;

D_val = [];
for i = 1:sv
    for j = 1:sh
        win_ss = binary_data_in( 1+(i-1)*win_sh:i*win_sh,1+(j-1)*win_sv:j*win_sv );
        counter_1 = numel(find(win_ss == 1));
        D_val(i,j) = counter_1/max_den;
    end
end
density_out = D_val(:);