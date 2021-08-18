function M0 = target_params(r_P_peaks,std_P_tg,P_peaks,range)
%
% Description:
%       Total power reflectivity (power received by a radar and integrated
%       over Doppler spectra) computed based on the target position, width
%       and maximum power, when consider that the target's power can be 
%       characterized by a Gaussian shape in the Doppler spectra.
%
% INPUT:
%       r_P_peaks       - vector, contains the distance [km] where the maximum 
%                         power of the target(s) is/are located;
%       std_P_tg        - vector, contains the width [km] of the target(s).
%                         This defines an extended/point target;
%       P_peaks         - vectors, contains the max power that each target
%                         has.
%       range           - vector, range bins.
%
% OUTPUT:
%       M0              - total power received for each range bin, when
%                         integrated over Doppler spectra.
%
% =======================================================================
% 18.08.2021 - AG, Helper added
%
% =======================================================================
K = numel(P_peaks);

% generate gaussian shape power distribution for each target
for i = 1:K    
    M0_tg(i,:) =  P_peaks(i)*exp(-(range-r_P_peaks(i)).^2/(2*std_P_tg(i)^2));
end
if K > 1
    M0 = sum(M0_tg);
else
    M0 = M0_tg;
end




