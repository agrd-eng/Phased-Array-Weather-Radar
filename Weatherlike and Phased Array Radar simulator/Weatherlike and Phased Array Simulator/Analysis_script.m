% Set the limits for FR main lobe and sidelobes
left_tg(1) = find(steer_theta<=8.799,1,'last')+1;
left_tg(2) = find(steer_theta>=10.77,1,'first');
right_tg(1) = find(steer_theta<=22.6,1,'last');
right_tg(2) = find(steer_theta>=26.5,1,'first');
extreme(1) = find(steer_theta<=2.88,1,'last')+1;
extreme(2) = find(steer_theta>=32.46,1,'first')-1;

xl_int = [steer_theta(left_tg(1)) steer_theta(left_tg(2))];
x2_int = [steer_theta(right_tg(1)) steer_theta(right_tg(2))];

xl_out1 = [steer_theta(extreme(1)) steer_theta(left_tg(1))];
xl_out2 = [steer_theta(left_tg(2)) steer_theta(right_tg(1))];
xl_out3 = [steer_theta(right_tg(2)) steer_theta(extreme(2))];

for j = 1:2 % 2 versions for Ground truth
    for i = 1:3% 3 moments
    M_truth(:,:,j) = [10*log10(abs(M0_truth(:,j)).');M1_truth(:,j).';M2_truth(:,j).'];
    end
end
M_est = zeros(3,numel(steer_theta),N_meth);

for m = 1:N_meth
    M_est(:,:,m) = [10*log10(abs(M0_est(m,:)));M1_est(m,:);M2_est(m,:)];
end


% M_rms_error(:,:,it) = squeeze(sqrt(var((M_est - repmat(M_truth,1,1,N_meth)),1,2)));
% M_rms_error(:,:,it) = sqrt(squeeze(sum(abs(M_est - repmat(M_truth,1,1,N_meth)).^2,2))./size(M_truth,2));



% -------------- ADD HERE ALL GRAPHS SCRIPTS YOU WANT TO PLOT -------------
analysis_graphs_1
% analysis_graphs_2_SNR
% -------------------------------------------------------------------------
