clc
clear
figure(1)
xText = 5.5;

load Run05_N10_NoiseFree
MSE_1PhaseWin9_NoiseFree(1,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(1,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N20_NoiseFree
MSE_1PhaseWin9_NoiseFree(2,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(2,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N30_NoiseFree
MSE_1PhaseWin9_NoiseFree(3,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(3,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N40_NoiseFree
MSE_1PhaseWin9_NoiseFree(4,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(4,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N50_NoiseFree
MSE_1PhaseWin9_NoiseFree(5,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(5,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N60_NoiseFree
MSE_1PhaseWin9_NoiseFree(6,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(6,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N70_NoiseFree
MSE_1PhaseWin9_NoiseFree(7,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(7,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N80_NoiseFree
MSE_1PhaseWin9_NoiseFree(8,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(8,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N90_NoiseFree
MSE_1PhaseWin9_NoiseFree(9,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(9,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N100_NoiseFree
MSE_1PhaseWin9_NoiseFree(10,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(10,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N110_NoiseFree
MSE_1PhaseWin9_NoiseFree(11,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(11,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N120_NoiseFree
MSE_1PhaseWin9_NoiseFree(12,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(12,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N130_NoiseFree
MSE_1PhaseWin9_NoiseFree(13,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(13,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N140_NoiseFree
MSE_1PhaseWin9_NoiseFree(14,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(14,1) = MSE_1phase_arccosinefree_withHarmonics(4);

load Run05_N150_NoiseFree
MSE_1PhaseWin9_NoiseFree(15,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9_NoiseFree(15,1) = MSE_1phase_arccosinefree_withHarmonics(4);

%% SNR=60dB
load Run05_N10
MSE_1PhaseWin9(1,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(1,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(1,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(1,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(1) = fs;


load Run05_N20
MSE_1PhaseWin9(2,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(2,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(2,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(2,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(2) = fs;

load Run05_N30
MSE_1PhaseWin9(3,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(3,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(3,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(3,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(3) = fs;

load Run05_N40
MSE_1PhaseWin9(4,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(4,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(4,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(4,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(4) = fs;

load Run05_N50
MSE_1PhaseWin9(5,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(5,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(5,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(5,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(5) = fs;

load Run05_N60
MSE_1PhaseWin9(6,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(6,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(6,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(6,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(6) = fs;

load Run05_N70
MSE_1PhaseWin9(7,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(7,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(7,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(7,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(7) = fs;

load Run05_N80
MSE_1PhaseWin9(8,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(8,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(8,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(8,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(8) = fs;

load Run05_N90
MSE_1PhaseWin9(9,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(9,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(9,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(9,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(9) = fs;

load Run05_N100
MSE_1PhaseWin9(10,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(10,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(10,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(10,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(10) = fs;

load Run05_N110
MSE_1PhaseWin9(11,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(11,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(11,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(11,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(11) = fs;

MSE_1PhaseWin9(12,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(12,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(12,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(12,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(12) = fs;

load Run05_N130
MSE_1PhaseWin9(13,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(13,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(13,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(13,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(13) = fs;

load Run05_N140
MSE_1PhaseWin9(14,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(14,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(14,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(14,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(14) = fs;

load Run05_N150
MSE_1PhaseWin9(15,1) = MSE_1phase_withHarmonics(4); % window length 9
MSE_1PhaseArccosinefreeWin9(15,1) = MSE_1phase_arccosinefree_withHarmonics(4);
MSE_1PhaseWin9_Filtered(15,1) = MSE_1phase_withHarmonics_filtered(4); 
MSE_1PhaseArccosinefreeWin9_Filtered(15,1) = MSE_1phase_arccosinefree_withHarmonics_filtered(4);
samplingFreq(15) = fs;

figure(1)
subplot(311)
plot(samplingFreq, MSE_1PhaseWin9_NoiseFree(1:15),'Color','black','LineWidth',1.5)
hold on
plot(samplingFreq , MSE_1PhaseArccosinefreeWin9_NoiseFree(1:15), '-.','Color','black','LineWidth',1.5)
ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
title('for noise-free signal','FontSize', fontSize, 'FontWeight', 'bold')
legend('1-phase', '1-phase arccosine-free')

subplot(312)
plot(samplingFreq, MSE_1PhaseWin9(1:15),'Color','black','LineWidth',1.5)
hold on
plot(samplingFreq , MSE_1PhaseArccosinefreeWin9(1:15), '-.','Color','black','LineWidth',1.5)
ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
title('for noisy signal, SNR = 60dB','FontSize', fontSize, 'FontWeight', 'bold')
legend('1-phase', '1-phase arccosine-free')

subplot(313)
plot(samplingFreq , MSE_1PhaseWin9_Filtered(1:15),'Color','black','LineWidth',1.5)
hold on
plot(samplingFreq , MSE_1PhaseArccosinefreeWin9_Filtered(1:15), '-.','Color','black','LineWidth',1.5)
xlabel('f_s','FontSize', fontSize, 'FontWeight', 'bold')
ylabel('MSE','FontSize', fontSize, 'FontWeight', 'bold')
title('for filtered noisy signal, SNR = 60dB','FontSize', fontSize, 'FontWeight', 'bold')
legend('1-phase', '1-phase arccosine-free')