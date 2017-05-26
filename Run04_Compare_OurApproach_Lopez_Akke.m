clc
clear
figure(1)
xText = 5.5;

load Run03Test_80dB
subplot(411)
bar(lengthWindow, MSE_1phase3phase_withHarmonics,'group')
ylabel('MSE','FontSize', fontSize, 'FontWeight', 'light')
yText = 0.9*max(max(MSE_1phase3phase_withHarmonics));
text(xText,yText, ['SNR = ' int2str(20*log10(1/sigma))])


load Run03Test_60dB
subplot(412)
bar(lengthWindow, MSE_1phase3phase_withHarmonics,'group')
ylabel('MSE','FontSize', fontSize, 'FontWeight', 'light')
yText = 0.9*max(max(MSE_1phase3phase_withHarmonics));
text(xText,yText, ['SNR = ' int2str(20*log10(1/sigma))])


load Run03Test_40dB
subplot(413)
bar(lengthWindow, MSE_1phase3phase_withHarmonics,'group')
ylabel('MSE','FontSize', fontSize, 'FontWeight', 'light')
yText = 0.9*max(max(MSE_1phase3phase_withHarmonics));
text(xText,yText, ['SNR = ' int2str(20*log10(1/sigma))])


load Run03Test_20dB
subplot(414)
bar(lengthWindow, MSE_1phase3phase_withHarmonics,'group')
ylabel('MSE','FontSize', fontSize, 'FontWeight', 'light')
yText = 0.9*max(max(MSE_1phase3phase_withHarmonics));
text(xText,yText, ['SNR = ' int2str(20*log10(1/sigma))])


gtext('[17]\downarrow')% Lopez
gtext('[17]\downarrow')% Lopez
gtext('[17]\downarrow')% Lopez
gtext('[17]\downarrow')% Lopez



gtext('[8]\downarrow') % Akke
gtext('[8]\downarrow') % Akke
gtext('[8]\downarrow') % Akke
gtext('[8]\downarrow') % Akke
