
clc
clear
close all

addpath ('D:\1_TUDelft\9_Master_Thesis\All programs\AG_finalver_software - Copy\Functions')
L = input('Give number of Doppler bins [256]: ');
if isempty(L)
    L = 256;
end
N = input('Give number of antenna elements [40]: ');
if isempty(N)
    N = 40;% number of antenna elements
end

% Fixed phased array parameter values
d = 20e-3; % distance between elements
f0 = 9650e6; %Hz center freq 
c = physconst('lightspeed');
lambda = c/f0; 
delta = d/lambda;% d/lambda (d = distance between antenna elements)
d_el = rad2deg(.886/(N*delta));% broadside elevation resolution
% d_el = 0.1;% if one wants to analyze Antenna Spatial Response
steer_theta = (-5:d_el:35);% BF elevation angles [deg]

% Max unambiguous Doppler velocity
vmax = 7.39;%m/s max unambigous Doppler velocity
dv = 2*vmax/L;
vel = -vmax:dv:vmax-dv;

% Target parameters
Prx_max = 10.^([40 60]./10);% represent power level at the mean theta position
el_Prx_max = [10 25]; % Elevation Angle that corresponds with targets max power
std_Prx_tg = 1.5*ones(1,numel(Prx_max));% Elevation Width of the spectrum

% % Generate M0
M0 = total_power_reflectivity(el_Prx_max,std_Prx_tg,Prx_max,steer_theta);

% Set noise variance w.r.t. the maximum power among all targets and given
% SNR
SNR = input('Give SNR value [dB]: [20]');
if isempty(SNR)
    SNR = 20;
end
varn = max(Prx_max)/10^(SNR/10);
comp = L*varn/N*dv;% compensation - computes noise floor for final plots

% M1 and M2
% Below is only one example to set these two parameters M1 and M2
m2 = 1; % Doppler spectrum width

M1 = linspace(-vmax+1*m2,vmax-1*m2,numel(M0)); % % monotonically increasing mean Doppler velocity
% or
% M1 = 5*ones(1,numel(M0)); % constant mean Doppler velocity
M2 = m2*ones(1,numel(M0)); % constant spectral width

% % Shows M0, M1 and M2
font_s = 30;
f1_1 = figure('units','normalized','outerposition',[0 0 1 1]);
figure(f1_1)

subplot(131)
plot(steer_theta, 10*log10(M0),'-ob','Linewidth',2)
ylim([-10 70])
xlim([0.5 29.5])
xlabel('Elevation [deg]')
ylabel('Power [dB]')
title('M0')
grid minor
set(gca,'Fontsize',font_s)

subplot(132)
plot(steer_theta, M1,'-or','Linewidth',2)
xlim([0.5 29.5])
xlabel('Elevation [deg]')
ylabel('Velocity [m/s]')
title('M1')
grid minor
set(gca,'Fontsize',font_s)

subplot(133)
plot(steer_theta, M2,'-og','Linewidth',2)
xlim([0.5 29.5])
xlabel('Elevation [deg]')
ylabel('Velocity [m/s]')
title('M2')
grid minor
set(gca,'Fontsize',font_s)

% Generate Weather like time domains signals to plug it in data model for
% the phased antenna array
[X,X_PS,M0_truth,M1_truth,M2_truth] = weatherlike_spectra_signals(M0,M1,M2,L,vmax,varn/N);
% Model used for data model of a phase array antenna is:
% X = S*X + Noise
% The unknown is X
S = gen_a(N,delta,steer_theta);
noise = sqrt(varn/2)*(randn(N,L)+1i*randn(N,L));% [antenna elem's]x[doppler bins]
Y = S*X+noise;
%
% Estimate X using three different beamforming methods:
% 1. recursive MMSE with gain control
% 2. Fourier, aka Matched Filtering (MF)
% 3. Capon, aka Minimum Variance Distortionless Response (MVDR)
N_meth = 3; % number of beamforming methods to compare
X_hat = zeros(numel(steer_theta),L,N_meth);% amplitude time domain signal
data_ED = zeros(numel(steer_theta),L,N_meth);% power spectrum signal

% Recursive MMSE (cMMSE), 'c' comes from converged, when th is reached
[X_hat(:,:,1),iter_out,W_cmmse] = cMMSE(Y,steer_theta,delta,varn,0.001);
% Fourier (FR)
[X_hat(:,:,2),W_fr] = FR(Y,steer_theta,delta);
% Capon (CP)
[X_hat(:,:,3),W_cp] = CP(Y,steer_theta,delta);


data_ED = 1/L.*abs(fftshift(fft(squeeze(X_hat(:,:,:)),[],2),2)).^2;% Power spectrum
for m = 1:N_meth
    
    for n = 1:numel(steer_theta)
        [M0_est(m,n),M1_est(m,n),M2_est(m,n)] = gauss_calc(data_ED(n,:,m),vmax);
    end

end
disp(['No of iterations for cMMSE: ',num2str(iter_out)])

% Code lines for plots are written in the next scripts
Analysis_script
test_spatial_response
estimation_accuracies

