
function M0 = total_power_reflectivity(el_Prx_max,std_Prx_tg,Prx_max,steer_theta)

%TOTAL_POWER_REFLECTIVITY Generates total power array for multiple targets
%
%   Description:
%   For fixed range and azimuth angle the atmospheric particles scanned in
%   one single radar volume is characterized by a total power value. This 
%   value can vary with elevation angle. This function returns vector M0
%   which is the superposition of the power reflected by multiple targets
%   at different elevation angles, discretized by the steering angle.
%   Target parameters are maximum power, position of this maximum w.r.t. to
%   elevation angle and the elevation width of the target's power
%   distribution. This function simulated targets power distribution vs
%   elevation angle as a Gaussian shape.
%
%   Usage:
%   M0 = total_power_reflectivity(el_Prx_max,std_Prx_tg,Prx_max,steer_theta)
%
%   Output:
%   M0          - vector, 1 x K, where K is the number of angles used for
%                 steering vector;
%   Input:
%   el_Prx_max  - vector, 1 x M, M is the number of targets 
%                 position of the maximum power w.r.t. elevation angle.
%   std_Prx_tg: - vector, 1 x M, 1st standard deviation of the target's power
%                 distribution; one can simulate either extended targets or
%                 point targets.               
%   Prx_max:    - vector, 1 x M, target's power (maximum value if we consider that each 
%                 target is characterized by a Gaussian distribution)
%

K = numel(Prx_max);

% generate gaussian shape power distribution for each target
for i = 1:K    
    M0_tg(i,:) =  Prx_max(i)*exp(-(steer_theta-el_Prx_max(i)).^2/(2*std_Prx_tg(i)^2));
end

if K > 1
    M0 = sum(M0_tg);
else
    M0 = M0_tg;
end

end

