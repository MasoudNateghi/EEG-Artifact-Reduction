%% load data
clear; close all; clc;
load("contaminated.mat")
load("pure.mat")
fs = 200; 
SrcName = cell(1, 19);
for i = 1:19
    SrcName{i} = strcat('s', num2str(i));
end
ChName = cell(1, 19);
for i = 1:19
    ChName{i} = strcat('ch', num2str(i));
end
%% part a
disp_eeg(contaminated, [], fs, ChName);
title('contaminated signal','Interpreter','Latex', 'FontSize', 10)
disp_eeg(pure, [], fs, ChName);
title('pure signal','Interpreter','Latex', 'FontSize', 10)
N = length(contaminated);
t = 1/fs:1/fs:N/fs;
T_on = find((t >= 1.35 & t <= 3.25) | (t >= 15.8 & t <= 16.57) | (t >= 19.47 & t <= 22.89));
t_on = (t >= 1.35 & t <= 3.25) | (t >= 15.8 & t <= 16.57) | (t >= 19.47 & t <= 22.89);
%% part b, c
C_tild_x = 1 / length(T_on) * (contaminated(:, T_on) * contaminated(:, T_on)');
C_x = 1 / 1e4 * (contaminated * contaminated');
[V, L] = eig(C_tild_x, C_x);
[~, index] = sort(diag(L), 'descend');
V = V(:, index);
S = V' * contaminated;
disp_eeg(S, [], fs, SrcName);
title('Sources for part c using GEVD','Interpreter','Latex', 'FontSize', 10)
SelSources = 3:19;
A = (V^(-1))';
X_den_GEVD = A(:, SelSources) * S(SelSources, :);
disp_eeg(X_den_GEVD, [], fs, ChName);
title('denoised signal using GEVD','Interpreter','Latex', 'FontSize', 10)
%% part d, e
% whitening the data
[coeff, score, latent] = pca((contaminated-mean(contaminated, 2))');
% coeff is the matrix of eigen vectors (sorted)
% score is zero mean whitened data
% latent has the eigen values
B = diag(latent)^(-1/2) * coeff';
z = diag(latent)^(-1/2) * score';    
d = 15;             % we use 15 principal components from 19
M = 2;              % the number of sources we find in DSS. here we extract two source
n = 100;            % 100 iteration for DSS termination
W = randn(d, M);    % Matrix of initial random w
for p = 1:M
    for j = 1:n
        wp = W(:, p);
        rp = wp' * z(1:d, :);
        rp_plus = rp .* t_on;
        wp_plus = sum(z(1:d, :) .* rp_plus, 2);
        if p ~= 1
            A = W(:, 1:p-1);
            wp_plus = (eye(d) - A*A') * wp_plus;
        end
        wp = wp_plus / norm(wp_plus, 2);
        W(:, p) = wp;
    end
end
r = W' * z(1:d, :);
disp_eeg(r, [], fs, []);
title('Sources for part e using DSS','Interpreter','Latex', 'FontSize', 10)
X_EOG = pinv(B(1:d, :)) * W * r + mean(contaminated, 2);
X_den_DSS = contaminated - X_EOG;
disp_eeg(X_den_DSS, [], fs, ChName);
title('denoised signal using DSS','Interpreter','Latex', 'FontSize', 10)
%% part f
RRMSE_GEVD = sqrt(sum(sum((pure - X_den_GEVD) .^ 2))) / sqrt(sum(sum(pure .^ 2)));
RRMSE_DSS  = sqrt(sum(sum((pure - X_den_DSS ) .^ 2))) / sqrt(sum(sum(pure .^ 2)));
RRMSE_worst_case = sqrt(sum(sum((pure - contaminated) .^ 2))) / sqrt(sum(sum(pure .^ 2)));






