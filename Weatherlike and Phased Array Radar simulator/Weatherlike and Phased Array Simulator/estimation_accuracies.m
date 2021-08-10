disp('--yellow--')
disp('---------------------')
disp('--total power--')
for i = 1:3
    temp = [10*log10(M0_truth(left_tg(1):left_tg(2),1).') - 10*log10(M0_est(i,left_tg(1):left_tg(2)))...
    10*log10(M0_truth(right_tg(1):right_tg(2),1).') - 10*log10(M0_est(i,right_tg(1):right_tg(2)))];
    mb_power(i) = mean(temp);
    std_power(i) = std(temp);
end
mb_power
std_power
disp('--mean Doppler velocity--')
for i = 1:3
    
    temp = [(M1_truth(left_tg(1):left_tg(2),1).') - (M1_est(i,left_tg(1):left_tg(2)))...
    (M1_truth(right_tg(1):right_tg(2),1).') - (M1_est(i,right_tg(1):right_tg(2)))];
    mb_vel(i) = mean(abs(temp));
    std_vel(i) = std(abs(temp));
end
mb_vel
std_vel
disp('--spectral width--')
for i = 1:3
    temp = [(M2_truth(left_tg(1):left_tg(2),1).') - (M2_est(i,left_tg(1):left_tg(2)))...
    (M2_truth(right_tg(1):right_tg(2),1).') - (M2_est(i,right_tg(1):right_tg(2)))];
    mb_sw(i) = mean(abs(temp));
    std_sw(i) = std(abs(temp));
end
mb_sw
std_sw

clear mb_power mb_vel mb_sw std_power std_vel std_sw

disp('----------------------------')
disp('---blue---')
disp('--total power--')
for i = 1:3
    temp = [10*log10(M0_truth(extreme(1):left_tg(1),1).') - 10*log10(M0_est(i,extreme(1):left_tg(1)))...
    10*log10(M0_truth(left_tg(2):right_tg(1),1).') - 10*log10(M0_est(i,left_tg(2):right_tg(1)))...
    10*log10(M0_truth(right_tg(2):extreme(2),1).') - 10*log10(M0_est(i,right_tg(2):extreme(2)))];
    mb_power(i) = mean(temp);
    std_power(i) = std(temp);
end
mb_power
std_power

disp('---mean D vel---')
for i = 1:3
    temp = [(M1_truth(extreme(1):left_tg(1),1).') - (M1_est(i,extreme(1):left_tg(1)))...
    (M1_truth(left_tg(2):right_tg(1),1).') - (M1_est(i,left_tg(2):right_tg(1)))...
    (M1_truth(right_tg(2):extreme(2),1).') - (M1_est(i,right_tg(2):extreme(2)))];
    mb_vel(i) = mean(abs(temp));
    std_vel(i) = std(abs(temp));
end
mb_vel
std_vel

disp('--spectral width---')
for i = 1:3
    temp = [(M2_truth(extreme(1):left_tg(1),1).') - (M2_est(i,extreme(1):left_tg(1)))...
    (M2_truth(left_tg(2):right_tg(1),1).') - (M2_est(i,left_tg(2):right_tg(1)))...
    (M2_truth(right_tg(2):extreme(2),1).') - (M2_est(i,right_tg(2):extreme(2)))];
    mb_sw(i) = mean(abs(temp));
    std_sw(i) = std(abs(temp));
end
mb_sw
std_sw
