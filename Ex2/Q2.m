%% Load Data
clear; close all; clc;
load("Ex2.mat")
load('Electrodes');
feq = 250 ;
ElecName = Electrodes.labels;
SrcName = cell(1, 32);
for i = 1:32
    SrcName{i} = strcat('s', num2str(i));
end
tempName = {'T3', 'FC5', 'T3', 'FC5', 'T3', 'FC5'};
%% X_noise_1 (SNR = -10dB)
SNR = -10;
Ps = sum(sum(X_org .^ 2));
Pn = sum(sum(X_noise_1 .^ 2));
sigma = sqrt(Ps / Pn * 10^(-SNR/10));

N = sigma*X_noise_1;
X = X_org + N;
offset = max(abs(X_org(:)));
disp_eeg(X, offset, feq, ElecName);
title('noisy EEG signal (X noise 1 SNR = -10dB)')
disp_eeg(X_org, offset, feq, ElecName);
title('original EEG signal')
%% PCA
N = size(X_org, 2);
mu_x = mean(X, 2);
Cx = (X - mu_x) * (X - mu_x)' / N;
[V, L] = eig(Cx);
[d, ind] = sort(diag(L), 'descend');
A = V(:, ind)';
B = diag(d) ^ (-1/2);
D = B*A;
y = D*X;
mu_y = mean(y, 2);
Cy = (y - mu_y) * (y - mu_y)' / N;
disp_eeg(y, offset, feq, SrcName);
title('sources of EEG signal (X noise 1 SNR = -10dB) using PCA')
%% Choose desired Components
y_den = zeros(32, 10240);
y_den(2, :) = y(2, :);
y_den(3, :) = y(3, :);
% y_den(4, :) = y(4, :);
y_den(5, :) = y(5, :);
y_den(7, :) = y(7, :);
y_den(9, :) = y(9, :);
y_den(16, :) = y(16, :);
% y_den(22, :) = y(22, :);
% y_den(23, :) = y(23, :);
y_den(32, :) = y(32, :);
X_den = pinv(D)*y_den;
disp_eeg(X_den, offset, feq, ElecName);
title('denoised EEG signal (X noise 1 SNR = -10dB) using PCA')
%% Plot Channel 13 & 24
temp = zeros(6, 10240);
temp(1:2, :) = X_org([13, 14], :);
temp(3:4, :) = X([13, 24], :);
temp(5:6, :) = X_den([13, 24], :);
disp_eeg(temp, offset, feq, tempName);
legend('original', 'original', 'noisy', 'noisy', 'denoised', 'denoised')
title('denoised channel 13 & 24 of EEG signal (X noise 1 SNR = -10dB) using PCA')
%% RRMSE
num = sum(sum((X_org - X_den) .^ 2));
denum = sum(sum(X_org .^ 2));
RRMSE = sqrt(num / denum)
num = sum(sum((X_org - X) .^ 2));
RRMSE = sqrt(num / denum)
%% ICA
[F,W,K]=COM2R(X,32);
y=pinv(F)*X;
disp_eeg(y, offset, feq, SrcName);
title('sources of EEG signal (X noise 1 SNR = -10dB) using ICA')
%% Choose desired Components
y_den = zeros(32, 10240);
y_den(3, :) = y(3, :);
y_den(20, :) = y(20, :);
y_den(25, :) = y(25, :);
y_den(31, :) = y(31, :);
X_den = F*y_den;
disp_eeg(X_den, offset, feq, ElecName);
title('denoised EEG signal (X noise 1 SNR = -10dB) using ICA')
%% Plot Channel 13 & 24
temp = zeros(6, 10240);
temp(1:2, :) = X_org([13, 14], :);
temp(3:4, :) = X([13, 24], :);
temp(5:6, :) = X_den([13, 24], :);
disp_eeg(temp, offset, feq, tempName);
legend('original', 'original', 'noisy', 'noisy', 'denoised', 'denoised')
title('denoised channel 13 & 24 of EEG signal (X noise 1 SNR = -10dB) using ICA')
%% RRMSE
num = sum(sum((X_org - X_den) .^ 2));
denum = sum(sum(X_org .^ 2));
RRMSE = sqrt(num / denum)
%% X_noise_1 (SNR = -20dB)
SNR = -20;
Ps = sum(sum(X_org .^ 2));
Pn = sum(sum(X_noise_1 .^ 2));
sigma = sqrt(Ps / Pn * 10^(-SNR/10));

N = sigma*X_noise_1;
X = X_org + N;
offset = max(abs(X_org(:)));
disp_eeg(X, offset, feq, ElecName);
title('noisy EEG signal (X noise 1 SNR = -20dB)')
disp_eeg(X_org, offset, feq, ElecName);
title('original EEG signal')
%% PCA
N = size(X_org, 2);
mu_x = mean(X, 2);
Cx = (X - mu_x) * (X - mu_x)' / N;
[V, L] = eig(Cx);
[d, ind] = sort(diag(L), 'descend');
A = V(:, ind)';
B = diag(d) ^ (-1/2);
D = B*A;
y = D*X;
mu_y = mean(y, 2);
Cy = (y - mu_y) * (y - mu_y)' / N;
disp_eeg(y, offset, feq, SrcName);
title('sources of EEG signal (X noise 1 SNR = -20dB) using PCA')
%% Choose desired Components
y_den = zeros(32, 10240);
y_den(2, :) = y(2, :);
y_den(3, :) = y(3, :);
y_den(4, :) = y(4, :);
y_den(5, :) = y(5, :);
y_den(7, :) = y(7, :);
y_den(9, :) = y(9, :);
y_den(16, :) = y(16, :);
y_den(22, :) = y(22, :);
y_den(23, :) = y(23, :);
y_den(32, :) = y(32, :);
X_den = pinv(D)*y_den;
disp_eeg(X_den, offset, feq, ElecName);
title('denoised EEG signal (X noise 1 SNR = -20dB) using PCA')
%% Plot Channel 13 & 24
temp = zeros(6, 10240);
temp(1:2, :) = X_org([13, 14], :);
temp(3:4, :) = X([13, 24], :);
temp(5:6, :) = X_den([13, 24], :);
offset = max(abs(temp(:)));
disp_eeg(temp, offset, feq, tempName);
legend('original', 'original', 'noisy', 'noisy', 'denoised', 'denoised')
title('denoised channel 13 & 24 of EEG signal (X noise 1 SNR = -20dB) using PCA')
%% RRMSE
num = sum(sum((X_org - X_den) .^ 2));
denum = sum(sum(X_org .^ 2));
RRMSE = sqrt(num / denum)
num = sum(sum((X_org - X) .^ 2));
RRMSE = sqrt(num / denum)
%% ICA
[F,W,K]=COM2R(X,32);
y=pinv(F)*X;
disp_eeg(y, offset, feq, SrcName);
title('sources of EEG signal (X noise 1 SNR = -20dB) using ICA')
%% Choose desired Components
y_den = zeros(32, 10240);
y_den(3, :) = y(3, :);
y_den(9, :) = y(9, :);
% y_den(11, :) = y(11, :);
y_den(22, :) = y(22, :);
y_den(24, :) = y(24, :);
y_den(30, :) = y(30, :);
y_den(31, :) = y(31, :);
X_den = F*y_den;
disp_eeg(X_den, offset, feq, ElecName);
title('denoised EEG signal (X noise 1 SNR = -20dB) using ICA')
%% Plot Channel 13 & 24
temp = zeros(6, 10240);
temp(1:2, :) = X_org([13, 14], :);
temp(3:4, :) = X([13, 24], :);
temp(5:6, :) = X_den([13, 24], :);
offset = max(abs(temp(:)));
disp_eeg(temp, offset, feq, tempName);
legend('original', 'original', 'noisy', 'noisy', 'denoised', 'denoised')
title('denoised channel 13 & 24 of EEG signal (X noise 1 SNR = -20dB) using ICA')
%% RRMSE
num = sum(sum((X_org - X_den) .^ 2));
denum = sum(sum(X_org .^ 2));
RRMSE = sqrt(num / denum)
%% X_noise_2 (SNR = -10dB)
SNR = -10;
Ps = sum(sum(X_org .^ 2));
Pn = sum(sum(X_noise_2 .^ 2));
sigma = sqrt(Ps / Pn * 10^(-SNR/10));

N = sigma*X_noise_2;
X = X_org + N;
offset = max(abs(X_org(:)));
disp_eeg(X, offset, feq, ElecName);
title('noisy EEG signal (X noise 2 SNR = -10dB)')
disp_eeg(X_org, offset, feq, ElecName);
title('original EEG signal')
%% PCA
N = size(X_org, 2);
mu_x = mean(X, 2);
Cx = (X - mu_x) * (X - mu_x)' / N;
[V, L] = eig(Cx);
[d, ind] = sort(diag(L), 'descend');
A = V(:, ind)';
B = diag(d) ^ (-1/2);
D = B*A;
y = D*X;
mu_y = mean(y, 2);
Cy = (y - mu_y) * (y - mu_y)' / N;
disp_eeg(y, offset, feq, SrcName);
title('sources of EEG signal (X noise 2 SNR = -10dB) using PCA')
%% Choose desired Components
y_den = zeros(32, 10240);
y_den(1, :) = y(1, :);
y_den(2, :) = y(2, :);
y_den(3, :) = y(3, :);
y_den(4, :) = y(4, :);
% y_den(15, :) = y(15, :);
% y_den(22, :) = y(22, :);
y_den(32, :) = y(32, :);
X_den = pinv(D)*y_den;
disp_eeg(X_den, offset, feq, ElecName);
title('denoised EEG signal (X noise 2 SNR = -10dB) using PCA')
%% Plot Channel 13 & 24
temp = zeros(6, 10240);
temp(1:2, :) = X_org([13, 14], :);
temp(3:4, :) = X([13, 24], :);
temp(5:6, :) = X_den([13, 24], :);
offset = max(abs(temp(:)));
disp_eeg(temp, offset, feq, tempName);
legend('original', 'original', 'noisy', 'noisy', 'denoised', 'denoised')
title('denoised channel 13 & 24 of EEG signal (X noise 2 SNR = -10dB) using PCA')
%% RRMSE
num = sum(sum((X_org - X_den) .^ 2));
denum = sum(sum(X_org .^ 2));
RRMSE = sqrt(num / denum)
num = sum(sum((X_org - X) .^ 2));
RRMSE = sqrt(num / denum)
%% ICA
[F,W,K]=COM2R(X,32);
y=pinv(F)*X;
disp_eeg(y, offset, feq, SrcName);
title('sources of EEG signal (X noise 2 SNR = -10dB) using ICA')
%% Choose desired Components
y_den = zeros(32, 10240);
y_den(6, :) = y(6, :);
y_den(25, :) = y(25, :);
X_den = F*y_den;
disp_eeg(X_den, offset, feq, ElecName);
title('denoised EEG signal (X noise 2 SNR = -10dB) using ICA')
%% Plot Channel 13 & 24
temp = zeros(6, 10240);
temp(1:2, :) = X_org([13, 14], :);
temp(3:4, :) = X([13, 24], :);
temp(5:6, :) = X_den([13, 24], :);
offset = max(abs(temp(:)));
disp_eeg(temp, offset, feq, tempName);
legend('original', 'original', 'noisy', 'noisy', 'denoised', 'denoised')
title('denoised channel 13 & 24 of EEG signal (X noise 2 SNR = -10dB) using ICA')
%% RRMSE
num = sum(sum((X_org - X_den) .^ 2));
denum = sum(sum(X_org .^ 2));
RRMSE = sqrt(num / denum)
%% X_noise_2 (SNR = -20dB)
SNR = -20;
Ps = sum(sum(X_org .^ 2));
Pn = sum(sum(X_noise_2 .^ 2));
sigma = sqrt(Ps / Pn * 10^(-SNR/10));

N = sigma*X_noise_2;
X = X_org + N;
offset = max(abs(X_org(:)));
disp_eeg(X, offset, feq, ElecName);
title('noisy EEG signal (X noise 2 SNR = -20dB)')
disp_eeg(X_org, offset, feq, ElecName);
title('original EEG signal')
%% PCA
N = size(X_org, 2);
mu_x = mean(X, 2);
Cx = (X - mu_x) * (X - mu_x)' / N;
[V, L] = eig(Cx);
[d, ind] = sort(diag(L), 'descend');
A = V(:, ind)';
B = diag(d) ^ (-1/2);
D = B*A;
y = D*X;
mu_y = mean(y, 2);
Cy = (y - mu_y) * (y - mu_y)' / N;
disp_eeg(y, offset, feq, SrcName);
title('sources of EEG signal (X noise 2 SNR = -20dB) using PCA')
%% Choose desired Components
y_den = zeros(32, 10240);
y_den(4, :) = y(4, :);
y_den(17, :) = y(17, :);
y_den(32, :) = y(32, :);
X_den = pinv(D)*y_den;
disp_eeg(X_den, offset, feq, ElecName);
title('denoised EEG signal (X noise 2 SNR = -20dB) using PCA')
%% Plot Channel 13 & 24
temp = zeros(6, 10240);
temp(1:2, :) = X_org([13, 14], :);
temp(3:4, :) = X([13, 24], :);
temp(5:6, :) = X_den([13, 24], :);
offset = max(abs(temp(:)));
disp_eeg(temp, offset, feq, tempName);
legend('original', 'original', 'noisy', 'noisy', 'denoised', 'denoised')
title('denoised channel 13 & 24 of EEG signal (X noise 2 SNR = -20dB) using PCA')
%% RRMSE
num = sum(sum((X_org - X_den) .^ 2));
denum = sum(sum(X_org .^ 2));
RRMSE = sqrt(num / denum)
num = sum(sum((X_org - X) .^ 2));
RRMSE = sqrt(num / denum)
%% ICA
[F,W,K]=COM2R(X,32);
y=pinv(F)*X;
offset = max(abs(y(:)));
disp_eeg(y, offset, feq, SrcName);
title('sources of EEG signal (X noise 2 SNR = -20dB) using ICA')
%% Choose desired Components
y_den = zeros(32, 10240);
y_den(21, :) = y(21, :);
X_den = F*y_den;
offset = max(abs(X_den(:)));
disp_eeg(X_den, offset, feq, ElecName);
title('denoised EEG signal (X noise 2 SNR = -20dB) using ICA')
%% Plot Channel 13 & 24
temp = zeros(6, 10240);
temp(1:2, :) = X_org([13, 14], :);
temp(3:4, :) = X([13, 24], :);
temp(5:6, :) = X_den([13, 24], :);
offset = max(abs(temp(:)));
disp_eeg(temp, offset, feq, tempName);
legend('original', 'original', 'noisy', 'noisy', 'denoised', 'denoised')
title('denoised channel 13 & 24 of EEG signal (X noise 2 SNR = -20dB) using ICA')
%% RRMSE
num = sum(sum((X_org - X_den) .^ 2));
denum = sum(sum(X_org .^ 2));
RRMSE = sqrt(num / denum)























