%% load dataset
clear; close all; clc;
load("Ex1_data.mat")
SrcName = cell(1, 8);
for i = 1:8
    SrcName{i} = strcat('s', num2str(i));
end
ChName = cell(1, 8);
for i = 1:8
    ChName{i} = strcat('ch', num2str(i));
end
fs = 100;
% visualze data
offset = max(abs(X_org(:)));
disp_eeg(X_org, offset, fs, ChName);
title('original signal', 'FontSize', 10)
t = 0:1/fs:(length(X_org)-1)/fs;
%% part a GEVD
P_x = 1 / (1e4 - 400) * (X_org(:, 1:end-400) * X_org(:, 401:end)'); 
C_x = 1 / 1e4 * (X_org * X_org');
P_tild_x = (P_x + P_x') / 2;
[V, ~] = eig(P_tild_x, C_x);
S = V' * X_org;
offset = max(abs(S(:)));
disp_eeg(S, offset, fs, SrcName);
title('Sources for part a using GEVD','Interpreter','Latex', 'FontSize', 10)
SelSources = 8;         % greatest eigen value is the last one. so we use 8th source
A = (V^(-1))';
X1_hat = A(:, SelSources) * S(SelSources, :);
offset = max(abs(X1_hat(:)));
disp_eeg(X1_hat, offset, fs, ChName);
title('$$\hat{x}_1$$(t) for part a using GEVD ','Interpreter','Latex', 'FontSize', 10)
offset = max(abs(X1(:)));
disp_eeg(X1, offset, fs, ChName);
title('$$x_1$$(t) for part a using GEVD ','Interpreter','Latex', 'FontSize', 10)
RRMSE = sqrt(sum(sum((X1 - X1_hat) .^ 2))) / sqrt(sum(sum(X1 .^ 2))); %#ok<NASGU>    
%% part b GEVD
interval = 300:700;
f = zeros(size(interval));
i = 1;
for tau = interval
    P_x = 1 / (1e4 - tau) * (X_org(:, 1:end-tau) * X_org(:, tau+1:end)'); 
    C_x = 1 / 1e4 * (X_org * X_org');
    P_tild_x = (P_x + P_x') / 2;
    [~, L] = eig(P_tild_x, C_x);
    f(i) = max(diag(L));
    i = i + 1;
end
[~, i] = max(f);
best_tau = interval(i);  
plot(interval, f)
xlabel('Time period')
ylabel('rayleigh quotient')

P_x = 1 / (1e4 - best_tau) * (X_org(:, 1:end-best_tau) * X_org(:, best_tau+1:end)'); 
C_x = 1 / 1e4 * (X_org * X_org');
P_tild_x = (P_x + P_x') / 2;
[V, ~] = eig(P_tild_x, C_x);
S = V' * X_org;
offset = max(abs(S(:)));
disp_eeg(S, offset, fs, SrcName);
title('Sources for part b using GEVD with best $$\tau$$ = 401', 'Interpreter','Latex', 'FontSize', 10)
SelSources = 8;         % greatest eigen value is the last one. so we use 8th source
A = (V^(-1))';
X1_hat = A(:, SelSources) * S(SelSources, :);
offset = max(abs(X1_hat(:)));
disp_eeg(X1_hat, offset, fs, ChName);
title('$$\hat{x}_1$$(t) for part b using GEVD with best $$\tau$$ = 401', 'Interpreter','Latex', 'FontSize', 10)
RRMSE = sqrt(sum(sum((X1 - X1_hat) .^ 2))) / sqrt(sum(sum(X1 .^ 2))); %#ok<NASGU>    
%% part c GEVD
t_on = find(T1 == 1);
C_tild_x = 1 / sum(T1) * (X_org(:, t_on) * X_org(:, t_on)');
C_x = 1 / 1e4 * (X_org * X_org');
[V, ~] = eig(C_tild_x, C_x);
S = V' * X_org;
offset = max(abs(S(:)));
disp_eeg(S, offset, fs, SrcName);
title('Sources for part c using GEVD','Interpreter','Latex', 'FontSize', 10)
SelSources = 8;
A = (V^(-1))';
X2_hat = A(:, SelSources) * S(SelSources, :);
offset = max(abs(X2_hat(:)));
disp_eeg(X2_hat, offset, fs, ChName);
title('$$\hat{x}_2$$(t) for part c using GEVD ','Interpreter','Latex', 'FontSize', 10)
offset = max(abs(X2(:)));
disp_eeg(X2, offset, fs, ChName);
title('$$x_2$$(t) for part c using GEVD ','Interpreter','Latex', 'FontSize', 10)
RRMSE = sqrt(sum(sum((X2 - X2_hat) .^ 2))) / sqrt(sum(sum(X2 .^ 2))); %#ok<NASGU> 
%% part d GEVD
t_on = find(T2 == 1);
C_tild_x = 1 / sum(T2) * (X_org(:, t_on) * X_org(:, t_on)');
C_x = 1 / 1e4 * (X_org * X_org');
[V, ~] = eig(C_tild_x, C_x);
S = V' * X_org;
figure;
plot(t, S(8, :))
hold on
plot(t, T2, "LineWidth", 2)
xlabel('t')
title('Sources for part d using GEVD with T2','Interpreter','Latex', 'FontSize', 10)
legend('source', 'T2')

T1_estimated = smoother(S(8, :), 2, 0.5, 2, 100);
figure;
plot(T1_estimated); hold on; plot(T1)
title('T1 vs. estimated T1')
legend('T1', 'T1 estimated')

t_on = find(T1_estimated == 1);
C_tild_x = 1 / sum(T1_estimated) * (X_org(:, t_on) * X_org(:, t_on)');
C_x = 1 / 1e4 * (X_org * X_org');
[V, ~] = eig(C_tild_x, C_x);
S = V' * X_org;
figure;
plot(t, S(8, :))
hold on
plot(t, T1_estimated, "LineWidth", 2)
xlabel('t')
title('Sources for part d using GEVD with estimated T1','Interpreter','Latex', 'FontSize', 10)
legend('source', 'estimated T1')

A = (V^(-1))';
SelSources = 8;
X2_hat = A(:, SelSources) * S(SelSources, :);
offset = max(abs(X2_hat(:)));
disp_eeg(X2_hat, offset, fs, ChName);
title('$$\hat{x}_2$$(t) for part d using GEVD with estimated T1','Interpreter','Latex', 'FontSize', 10)
RRMSE = sqrt(sum(sum((X2 - X2_hat) .^ 2))) / sqrt(sum(sum(X2 .^ 2))); %#ok<NASGU> 
%% part e GEVD
X = fft(X_org')';
X = [X(:, 1e4/2+1:end) X(:, 1:1e4/2)];
N = length(X_org);
k = 0:floor((N - 1) / 2);
fk_positive = k * fs / N;
k = floor((N - 1) / 2):N-1;
fk_negative = -(N - k) / N * fs;
f = [fk_negative(2:end) fk_positive];
nu = find(f >= -15 & f <=-10 | f >= 10 & f <= 15);
S_x = 1 / length(nu) * (X(:, nu)*X(:, nu)');
C_x = 1 / 1e4 * (X_org * X_org'); 
[V, ~] = eig(S_x, C_x);
V = real(V);
Sf = V' * X;

figure;
for i = 1:8
    subplot(4, 2, i)
    plot(f, abs(Sf(i, :)))
    xlabel('Frequency (Hz)')
    str = strcat('Source ', num2str(i));
    ylabel(str)
end
sgtitle('Sources for part e using GEVD ','Interpreter','Latex', 'FontSize', 10)

S = V' * X_org;
SelSources = 1;
A = (V^(-1))';
X3_hat = A(:, SelSources) * S(SelSources, :);
X3_hat_f = fft(X3_hat')';
X3_hat_f = [X3_hat_f(:, 1e4/2+1:end) X3_hat_f(:, 1:1e4/2)];

figure;
for i = 1:8
    subplot(4, 2, i)
    plot(f, abs(X3_hat_f(i, :)))
    xlabel('Frequency (Hz)')
    str = strcat('Source ', num2str(i));
    ylabel(str)
end
sgtitle('$$\hat{X}_3$$(f) for part e using GEVD ','Interpreter','Latex', 'FontSize', 10)

X3_f = fft(X3')';
X3_f = [X3_f(:, 1e4/2+1:end) X3_f(:, 1:1e4/2)];
figure;
for i = 1:8
    subplot(4, 2, i)
    plot(f, abs(X3_f(i, :)))
    xlabel('Frequency (Hz)')
    str = strcat('Source ', num2str(i));
    ylabel(str)
end
sgtitle('$$X_3$$(f) for part e using GEVD ','Interpreter','Latex', 'FontSize', 10)

RRMSE = sqrt(sum(sum((X3 - X3_hat) .^ 2))) / sqrt(sum(sum(X3 .^ 2))); %#ok<NASGU> 
%% part f GEVD
nu = find(f >= -25 & f <=-5 | f >= 5 & f <= 25);
S_x = 1 / length(nu) * (X(:, nu)*X(:, nu)');
C_x = 1 / 1e4 * (X_org * X_org'); 
[V, ~] = eig(S_x, C_x);
V = real(V);
Sf = V' * X;
f_on2 = f >= -25 & f <=-5 | f >= 5 & f <= 25;
f_on2 = double(f_on2);

figure;
offset = max(abs(Sf(1, :)));
plot(f, abs(Sf(1, :))); hold on; plot(f, offset*f_on2);
xlabel('Frequency (Hz)')
ylabel('Source 1')
title('Sources for part f using GEVD with f = [5 25]','Interpreter','Latex', 'FontSize', 10)
legend('source', 'active frequencies')

f_on_estimated = smoother(abs(Sf(1, :)), 100, 0.5, 10, fs);
nu = f_on_estimated == 1;
S_x = 1 / length(nu) * (X(:, nu)*X(:, nu)');
C_x = 1 / 1e4 * (X_org * X_org'); 
[V, ~] = eig(S_x, C_x);
V = real(V);
Sf = V' * X;

figure;
offset = max(abs(Sf(1, :)));
plot(f, abs(Sf(1, :))); hold on; plot(f, offset*f_on_estimated);
xlabel('Frequency (Hz)')
ylabel('Source 1')
title('Sources for part f using GEVD with with estimated f','Interpreter','Latex', 'FontSize', 10)
legend('source', 'estimated active frequencies')

S = V' * X_org;
SelSources = 1;
A = (V^(-1))';
X3_hat = A(:, SelSources) * S(SelSources, :);
X3_hat_f = fft(X3_hat')';
X3_hat_f = [X3_hat_f(:, 1e4/2+1:end) X3_hat_f(:, 1:1e4/2)];

figure;
for i = 1:8
    subplot(4, 2, i)
    plot(f, abs(X3_hat_f(i, :)))
    xlabel('Frequency (Hz)')
    str = strcat('Source ', num2str(i));
    ylabel(str)
end
sgtitle('$$\hat{X}_3$$(f) for part f using GEVD with estimated f','Interpreter','Latex', 'FontSize', 10)

RRMSE = sqrt(sum(sum((X3 - X3_hat) .^ 2))) / sqrt(sum(sum(X3 .^ 2))); %#ok<NASGU> 
%% whitening the data
[coeff, score, latent] = pca((X_org-mean(X_org, 2))');
% coeff is the matrix of eigen vectors (sorted)
% score is zero mean whitened data
% latent has the eigen values
B = diag(latent)^(-1/2) * coeff';
z = diag(latent)^(-1/2) * score';    
%% part a DSS
d = 6;      % we use six principal components from eight
M = 1;      % the number of sources we find in DSS. here we extract one source
n = 1000;   % 1000 iteration for DSS termination
W = randn(d, M);    % Matrix of initial random w
for p = 1:M
    for j = 1:n
        wp = W(:, p);           % initialization for wp          
        rp = wp' * z(1:d, :);   % noisy estimation of the source
        L = 1e4 / 400;          % number of periods
        r_tild = zeros(1, 400);
        for l = 1:L
            r_tild = r_tild + rp(1+(l-1)*400:400*l) / L;    % average over all periods
        end
        rp_plus = repmat(r_tild, [1 L]);                    % f(rp) 
        wp_plus = sum(z(1:d, :) .* rp_plus, 2);             % ML estimation
        if p ~= 1                                           % orthogonalization
            A = W(:, 1:p-1);
            wp_plus = (eye(d) - A*A') * wp_plus;
        end
        wp = wp_plus / norm(wp_plus, 2);                    % normalization
        W(:, p) = wp;                                       % update 
    end
end
r = W' * z(1:d, :);         % we do not extract all sources so we use this line of code
plot(t, r)
xlabel('time (seconds)')
ylabel('source 1')
title('Sources for part a using DSS','Interpreter','Latex', 'FontSize', 10)
X1_hat = pinv(B(1:d, :)) * W * r + mean(X_org, 2);  % display x1_hat
offset = max(abs(X1_hat(:)));
disp_eeg(X1_hat, offset, fs, ChName);
title('$$\hat{x}_1$$(t) for part a using DSS ','Interpreter','Latex', 'FontSize', 10)
RRMSE = sqrt(sum(sum((X1 - X1_hat) .^ 2))) / sqrt(sum(sum(X1 .^ 2))); %#ok<NASGU> 
%% part b DSS
%  !!!!!!!!!! LONG RUN TIME !!!!!!!!!!!!!
RRMSE = zeros(size(interval));
i = 1; n = 100;
for tau = interval
    W = randn(d, M);
    for p = 1:M
        for j = 1:n
            wp = W(:, p);
            rp = wp' * z(1:d, :);
            L = floor(1e4 / tau);
            r_tild = zeros(1, tau);
            for l = 1:L
            r_tild = r_tild + rp(1+(l-1)*tau:tau*l) / L;
                if l == L
                    tmp1 = rp(tau*l+1: end);
                    tmp3 = tau - rem(1e4, tau);
                    tmp2 = r_tild(end-tmp3+1:end);
                    tmp = [tmp1, tmp2];
                    r_tild = r_tild + tmp / L;
                end
            end
            main_part = repmat(r_tild, [1 L]);
            remaining_part = r_tild(1:rem(1e4, tau));
            rp_plus = [main_part, remaining_part];
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
    X1_hat = pinv(B(1:d, :)) * W * r + mean(X_org, 2);
    RRMSE(i) = sqrt(sum(sum((X1 - X1_hat) .^ 2))) / sqrt(sum(sum(X1 .^ 2)));
    i = i + 1;
    tau %#ok<NOPTS> 
end

figure;
plot(interval, RRMSE)
xlabel('period'); ylabel('RRMSE')
[~, i] = min(RRMSE);
best_tau = interval(i);

xlabel('Time period')
ylabel('rayleigh quotient')
W = randn(d, M);
for p = 1:M
    for j = 1:n
        wp = W(:, p);
        rp = wp' * z(1:d, :);
        L = floor(1e4 / best_tau);
        r_tild = zeros(1, best_tau);
        for l = 1:L
        r_tild = r_tild + rp(1+(l-1)*best_tau:best_tau*l) / L;
            if l == L
                tmp1 = rp(best_tau*l+1: end);
                tmp3 = best_tau - rem(1e4, best_tau);
                tmp2 = r_tild(end-tmp3+1:end);
                tmp = [tmp1, tmp2];
                r_tild = r_tild + tmp / L;
            end
        end
        main_part = repmat(r_tild, [1 L]);
        remaining_part = r_tild(1:rem(1e4, best_tau));
        rp_plus = [main_part, remaining_part];
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
figure;
plot(t, r)
xlabel('time (seconds)')
ylabel('source 1')
title('Sources for part b using DSS with best $$\tau$$ = 395','Interpreter','Latex', 'FontSize', 10)
X1_hat = pinv(B(1:d, :)) * W * r + mean(X_org, 2);  
offset = max(abs(X1_hat(:)));
disp_eeg(X1_hat, offset, fs, ChName);
title('$$\hat{x}_1$$(t) for part b using DSS with best $$\tau$$ = 395','Interpreter','Latex', 'FontSize', 10)
RRMSE = sqrt(sum(sum((X1 - X1_hat) .^ 2))) / sqrt(sum(sum(X1 .^ 2))); %#ok<NASGU> 
%% part c DSS
d = 6;
M = 1; n = 1000; 
z = diag(latent)^(-1/2) * score';
W = randn(d, M);
for p = 1:M
    for j = 1:n
        wp = W(:, p);
        rp = wp' * z(1:d, :);
        rp_plus = rp .* T1;
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
plot(t, r); hold on; plot(t, T1, 'LineWidth', 2);
title('Sources for part c using DSS','Interpreter','Latex', 'FontSize', 10)
X2_hat = pinv(B(1:d, :)) * W * r + mean(X_org, 2);
offset = max(abs(X2_hat(:)));
disp_eeg(X2_hat, offset, fs, ChName);
title('$$\hat{x}_2$$(t) for part c using DSS ','Interpreter','Latex', 'FontSize', 10)
RRMSE = sqrt(sum(sum((X2 - X2_hat) .^ 2))) / sqrt(sum(sum(X2 .^ 2))); %#ok<NASGU> 
%% part d DSS
d = 6;
M = 1; n = 1000; 
W = randn(d, M);
for p = 1:M
    for j = 1:n
        wp = W(:, p);
        rp = wp' * z(1:d, :);
        rp_plus = rp .* T2;
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
figure;
plot(t, r)
hold on
plot(t, T2, "LineWidth", 2)
xlabel('t')
title('Sources for part d using DSS with T2','Interpreter','Latex', 'FontSize', 10)
%
T1_estimated = smoother(r, 2, 0.5, 2, 100);
figure;
plot(T1_estimated); hold on; plot(T1)
title('T1 vs. estimated T1')
legend('T1', 'T1 estimated')
%
d = 6;
M = 1; n = 1000; 
W = randn(d, M);
for p = 1:M
    for j = 1:n
        wp = W(:, p);
        rp = wp' * z(1:d, :);
        rp_plus = rp .* T1_estimated;
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
figure;
plot(t, r)
hold on
plot(t, T1_estimated, "LineWidth", 2)
xlabel('t')
title('Sources for part d using DSS with estimated T1','Interpreter','Latex', 'FontSize', 10)
X2_hat = pinv(B(1:d, :)) * W * r + mean(X_org, 2);
offset = max(abs(X2_hat(:)));
disp_eeg(X2_hat, offset, fs, ChName);
title('$$\hat{x}_2$$(t) for part c using DSS with estimated T1','Interpreter','Latex', 'FontSize', 10)
RRMSE = sqrt(sum(sum((X2 - X2_hat) .^ 2))) / sqrt(sum(sum(X2 .^ 2))); %#ok<NASGU> 
%% part e DSS
d = 6;
M = 1; n = 100; 
z = diag(latent)^(-1/2) * score';
W = randn(d, M);
for p = 1:M
    for j = 1:n
        wp = W(:, p);
        rp = wp' * z(1:d, :);
        rp_plus = bandpass(rp, [10 15], fs);
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
r_f = fft(r);
r_f = [r_f(:, 1e4/2+1:end) r_f(:, 1:1e4/2)];
figure;
plot(f, abs(r_f))
title('Sources for part e using DSS','Interpreter','Latex', 'FontSize', 10)
X3_hat = pinv(B(1:d, :)) * W * r + mean(X_org, 2);
X3_hat_f = fft(X3_hat')';
X3_hat_f = [X3_hat_f(:, 1e4/2+1:end) X3_hat_f(:, 1:1e4/2)];
figure;
for i = 1:8
    subplot(4, 2, i)
    plot(f, abs(X3_hat_f(i, :)))
    xlabel('Frequency (Hz)')
    str = strcat('Source ', num2str(i));
    ylabel(str)
end
sgtitle('$$\hat{X}_3$$(f) for part e using DSS ','Interpreter','Latex', 'FontSize', 10)
RRMSE = sqrt(sum(sum((X3 - X3_hat) .^ 2))) / sqrt(sum(sum(X3 .^ 2))); %#ok<NASGU> 
%% part f DSS
d = 6;
M = 1; n = 100; 
z = diag(latent)^(-1/2) * score';
W = randn(d, M);
for p = 1:M
    for j = 1:n
        wp = W(:, p);
        rp = wp' * z(1:d, :);
        rp_plus = bandpass(rp, [5 25], fs);
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
r_f = fft(r);
r_f = [r_f(:, 1e4/2+1:end) r_f(:, 1:1e4/2)];
figure;
offset = max(abs(r_f(1, :)));
plot(f, abs(r_f)); hold on; plot(f, offset*f_on2);
xlabel('Frequency (Hz)')
ylabel('Source 1')
title('Sources for part f using DSS with f = [5 25]','Interpreter','Latex', 'FontSize', 10)
legend('source', 'active frequencies')
f_on_estimated = smoother(abs(r_f(1, :)), 100, 0.4, 15, fs);
% extract positive frequency interval
tmp = f(f_on_estimated > 0);
f_min = min(tmp(tmp > 0));
f_max = max(tmp(tmp > 0));
d = 6;
M = 1; n = 100; 
z = diag(latent)^(-1/2) * score';
W = randn(d, M);
for p = 1:M
    for j = 1:n
        wp = W(:, p);
        rp = wp' * z(1:d, :);
        rp_plus = bandpass(rp, [f_min f_max], fs);
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
r_f = fft(r);
r_f = [r_f(:, 1e4/2+1:end) r_f(:, 1:1e4/2)];
figure;
offset = max(abs(r_f(1, :)));
plot(f, abs(r_f)); hold on; plot(f, offset*f_on_estimated);
xlabel('Frequency (Hz)')
ylabel('Source 1')
title('Sources for part f using DSS with estimated f','Interpreter','Latex', 'FontSize', 10)
legend('source', 'active frequencies')
X3_hat = pinv(B(1:d, :)) * W * r + mean(X_org, 2);
X3_hat_f = fft(X3_hat')';
X3_hat_f = [X3_hat_f(:, 1e4/2+1:end) X3_hat_f(:, 1:1e4/2)];
figure;
for i = 1:8
    subplot(4, 2, i)
    plot(f, abs(X3_hat_f(i, :)))
    xlabel('Frequency (Hz)')
    str = strcat('Source ', num2str(i));
    ylabel(str)
end
sgtitle('$$\hat{X}_3$$(f) for part f using DSS ','Interpreter','Latex', 'FontSize', 10)
RRMSE = sqrt(sum(sum((X3 - X3_hat) .^ 2))) / sqrt(sum(sum(X3 .^ 2))); 













