%% Load Data
clear; close all; clc; 
load("NewData3.mat")
load("Electrodes.mat")
SrcName = cell(1, 32);
for i = 1:32
    SrcName{i} = strcat('s', num2str(i));
end
feq = 250; % sampling frequency
elabels = Electrodes.labels;
elocsX = Electrodes.X;
elocsY = Electrodes.Y;
%% Visualize Data
offset = max(abs(EEG_Sig(:))); % NewData1, 2, 3
% offset = max(abs(EEG_Sig(:)))/100; %NewData 4
disp_eeg(EEG_Sig, offset, feq, elabels);
title('original EEG signal')
%% ICA
[F,W,K]=COM2R(EEG_Sig,21);
y=pinv(F)*EEG_Sig;
offset = max(abs(y(:)));
disp_eeg(y, offset, feq, SrcName);
title('sorces of EEG signal')
%% 
N = size(EEG_Sig, 2);
t = 0:1/feq:(N-1)/feq;
for i = 1:21
    figure;
    subplot(3, 1, 1)
    plot(t, y(i, :))
    title(strcat('source ', num2str(i)))
    subplot(3, 1, 2)
    L = 512; 
    noverlap = 256; 
    nfft = L;
    [pxx, f] = pwelch(y(i, :), L, noverlap, nfft, feq);
    plot(f, 10*log10(pxx))
    subplot(3, 1, 3)
    plottopomap(elocsX,elocsY,elabels,F(:, i))
end
%% Select Desired Sources 
% SelSources = [1:3, 5:9, 11:21]; % NewData1.mat
% SelSources = [1:5, 8:21]; % NewData2.mat
SelSources = [2:6, 8:21]; % NewData3.mat
% SelSources = [9, 20, 21]; % NewData4.mat
X_den = F(:, SelSources) * y(SelSources, :);
offset = max(abs(X_den(:)));
disp_eeg(X_den, offset, feq, elabels);
title('denoised EEG signal')





