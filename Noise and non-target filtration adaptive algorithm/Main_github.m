addpath('D:\1_TUDelft\9_Master_Thesis\All programs\AG_finalver_software - Copy\To_commit_2\Functions')
load('noise-range-average-power.mat')
% A)

% % Rate of false alarme as function of threshold [db]
Th = 0; % dB threshold used for noise clipping
win_size = 2; % window size used for AMF algorithm
rmse_th = 1e-3; % threshold used for AMF

% % SNR, range bin of the target and target width
% NOTE: number of targets can be set by giving multiple values for the
% following three variables
SNR_db = 30;
r_P_peaks = 7; % range [km] correspond with max Power of each cluster 
std_M0 = 0.5; % km

% % Radar Parameters
% Doppler bins and velocity axis
L = 256; 
vmax = 7.39;%m/s max unambigous Doppler velocity
dv = 2*vmax/L;
v_axis = -vmax:dv:vmax-dv;

% Range bins
Nr = numel(noise_range_average);
c = physconst('lightspeed');
B = 66e6;      %Hz
dr = c/(2*B)*1e-3;
range = (1:Nr)*dr;%[km]

% % Compute moments
% M0 - total power reflectivity
SNR_peaks = 10.^(SNR_db./10);% linear scale
Noise_var_peaks  = noise_range_average(round(r_P_peaks./dr+1));
P_peaks = SNR_peaks.*Noise_var_peaks;
std_P_tg = std_M0*ones(1,numel(r_P_peaks));% Range Width of each Cloud

M0 = target_params(r_P_peaks,std_P_tg,P_peaks,range);
% M1 and M2
m2 = .5; % Doppler spectrum width
M1 = linspace(-vmax+1*m2,vmax-1*m2,Nr); % % monotonically increasing mean Doppler velocity
% M1 = 5*ones(1,numel(M0)); % constant mean Doppler velocity
M2 = m2*ones(1,numel(M0)); % spectral width


% % Weatherlike time domain signals and ground truth moments
[X_PS_nf, X_PS_n,M0_truth,M1_truth,M2_truth] = weather_sig_simulator_beta(M0,M1,M2,L,vmax,noise_range_average);
% -------------------------------------------------------------------------

% B) Create ground truth
TG_clip_NF = zeros(size(X_PS_nf));
noise_var_mat = repmat(noise_range_average,[1, L]);
TG_clip_NF(abs(X_PS_nf) > noise_var_mat) = 1;

E_se = strel('square',2);
TG_clip_MF = imdilate(TG_clip_NF,E_se);

X_PS_nf_C = TG_clip_MF.*X_PS_nf;
% -------------------------------------------------------------------------

% C) AMF
% CLIPPING
TG_clip = zeros(size(X_PS_n));
noise_var_mat = repmat(noise_range_average*(10^(Th/10)),[1, L]);
TG_clip(abs(X_PS_n) > noise_var_mat) = 1;% Raw binary image
X_PS_n_RB = TG_clip .* X_PS_n; % image after noise thresholding

rmse_diff = []; % contains square diference between consecutive rmse values 
final_dim_SE = []; % final dimension of the structural element (SE)
 
av_ov = 0;% cst used to avoid overshooting
mem = 0;% cst to save the number of consecutive H1 hypothesis
rep = 3;% threshold for the maximum number of consecutive H1
 
k = 1;
dim_SE = 1;
TG_clip_MF = []; % Target Mask
 
while 1
E_se = strel('square',dim_SE);
TG_clip_MF_temp = imopen(TG_clip,E_se);
TG_clip_MF(:,:,k) = imdilate(TG_clip_MF_temp,E_se);
% figure(1)
% imagesc(squeeze(TG_clip_MF(:,:,k)));
density = binary_scan_win(squeeze(TG_clip_MF(:,:,k)),win_size);
rmse(k) = density_pdf_fitting(density);

% --------------------AMF-----------------
if k > 1
       rmse_diff = [rmse_diff (abs(rmse(k)-rmse(k-1)))^2];
    if numel(final_dim_SE) == 0
        if (rmse(k)-rmse(k-1))^2 <= rmse_th
            if mem == k-1 || mem == 0                
                av_ov = av_ov+1;
                mem = k;
            else
                av_ov = 0;
                mem = 0;
            end
            if av_ov > rep
                break
            end
        else
            mem = 0;
        end   
    end
end

k = k+1;
dim_SE = dim_SE+1;
end
dim_SE_final = k-(rep);
disp(['Final SE dim: ',num2str(dim_SE_final)])

X_PS_n_C = squeeze(TG_clip_MF(:,:,dim_SE_final)).*X_PS_n;
% ------------------------------------------------------------------------

% D Plot

% ---
figure
subplot(221)
imagesc(v_axis,range,10*log10(abs(X_PS_n)))
title('Raw image')
xlabel('Doppler velocity [m/s]')
ylabel('Range [km]')
set(gca,'ydir','norm')

subplot(222)
imagesc(v_axis,range,10*log10(abs(X_PS_nf_C)))
title('Ground truth image')
xlabel('Doppler velocity [m/s]')
ylabel('Range [km]')
set(gca,'ydir','norm')

subplot(223)
imagesc(v_axis,range,10*log10(abs(X_PS_n_RB)))
title({'Image after noise thresholding',['Threshold: ', num2str(Th),' dB']})
xlabel('Doppler velocity [m/s]')
ylabel('Range [km]')
set(gca,'ydir','norm')

subplot(224)
imagesc(v_axis,range,10*log10(abs(X_PS_n_C)))
title('Image after AMF')
xlabel('Doppler velocity [m/s]')
ylabel('Range [km]')
set(gca,'ydir','norm')

%----
figure
plot([1:numel(rmse_diff)],rmse_diff,'-*r'),hold on
yline(rmse_th,'k','linewidth',2);hold off
ylim([rmse_th-0.5 +inf ])
title('AMF Converging curve')
ylabel('(rmse(k)-rmse(k-1))^2')
xlabel('SE dimension')

%-----
% Between noisy and ground truth image
for m = 1:Nr
    [M0_truth(m,2),M1_truth(m,2),M2_truth(m,2)] = gauss_calc(abs(X_PS_nf_C(m,:)),vmax,L);% Moments computation with noise
end
M0_truth(M0_truth == 0) = NaN;
% 1 - noisy
% 2 - noise free
M0_e(:,1) = sqrt((M0_truth(:,1) - M0_truth(:,2)).^2);
M1_e(:,1) = sqrt((M1_truth(:,1) - M1_truth(:,2)).^2);
M2_e(:,1) = sqrt((M2_truth(:,1) - M2_truth(:,2)).^2);

% Between noisy and after-AMF-image
for m = 1:Nr
    [M0_truth(m,3),M1_truth(m,3),M2_truth(m,3)] = gauss_calc(abs(X_PS_n_C(m,:)),vmax,L);% Moments computation with noise
end
M0_truth(M0_truth == 0) = NaN;

M0_e(:,2) = sqrt((M0_truth(:,3) - M0_truth(:,2)).^2);
M1_e(:,2) = sqrt((M1_truth(:,3) - M1_truth(:,2)).^2);
M2_e(:,2) = sqrt((M2_truth(:,3) - M2_truth(:,2)).^2);

M0_err_cor = nanmean(log10(M0_e(:,1)./M0_e(:,2)));
M1_err_cor = nanmean(log10(M1_e(:,1)./M1_e(:,2)));
M2_err_cor = nanmean(log10(M2_e(:,1)./M2_e(:,2)));

disp('-------------------------------------------')
disp('Estimation Accuracy Gain is: ')
disp(['M0: ',num2str(M0_err_cor),' [log10 scale]'])
disp(['M1: ',num2str(M1_err_cor),' [log10 scale]'])
disp(['M2: ',num2str(M2_err_cor),' [log10 scale]'])
