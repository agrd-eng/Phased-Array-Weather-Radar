% Test spatial response

% Gives results as in article when next variables are set in Main as follows:
% N = 128
% d_el = 0.1;
% steer_theta = (10:d_el:35);% BF elevation angles [deg]
% Prx_max = 10.^([20 40]./10);% represent power level at the mean theta position
% el_Prx_max = [25 60]; % Elevation Angle that corresponds with targets max power
% std_Prx_tg = 1*ones(1,numel(Prx_max));% Elevation Width of the spectrum
% SNR = 50

A = exp((0:N-1)'*1i*(2*pi*delta*sind(steer_theta)));
% ind(1)=find(steer_theta<el_Prx_max(1),1,'last')+1;
% ind(2)=find(steer_theta<el_Prx_max(2),1,'last')+1;
ind(1)=find(steer_theta<5,1,'last')+1;
% ind(2)=find(steer_theta<el_Prx_max(2),1,'last')+1;
% FR
y_FR = W_fr(:,ind(1))'*A;
% Cmmse
y_cmmse = W_cmmse(:,ind(1))'*A;
% CP
y_cp = W_cp(:,ind(1))'*A;

figure
plot(steer_theta,20*log10(abs(y_FR)),'b','linewidth',2),hold on
plot(steer_theta,20*log10(abs(y_cmmse)),'r','linewidth',2)
plot(steer_theta,20*log10(abs(y_cp)),'g','linewidth',2)
xline(steer_theta(ind(1)),':','linewidth',2);
% xline(steer_theta(ind(2)),':','linewidth',2);
title('Antenna Spatial Response')
xlabel('Elevation [deg]')
ylabel('Normalized directionality [dB]')

legend('FR','cMMSE','CP')
grid minor
% ylim([-70 40])
set(gca,'fontsize',20)
