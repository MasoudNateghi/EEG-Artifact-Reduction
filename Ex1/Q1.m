%% Load Data
clear; close all; clc;
load("Ex1.mat");
%% Visualize Data
figure;
scatter3(X(1, :), X(2, :), X(3, :))
xlim([-2 2])
ylim([-40 40])
zlim([-0.4 0.4])
xlabel('x1')
ylabel('x2')
zlabel('x3')
figure;
scatter3(X(1, :), X(2, :), X(3, :))
xlim([-40 40])
ylim([-40 40])
zlim([-40 40])
xlabel('x1')
ylabel('x2')
zlabel('x3')
title('original data')
%% Skewness
% % first feature
% m1 = mean(X(1, :));
% m2 = mean(X(1, :) .^ 2);
% m3 = mean(X(1, :) .^ 3);
% k2 = m2 - m1^2;
% k3 = m3 - 3*m2*m1 + 2*m1^3;
% k3_tild = k3 / (k2)^(3/2)
% % second feature
% m1 = mean(X(2, :));
% m2 = mean(X(2, :) .^ 2);
% m3 = mean(X(2, :) .^ 3);
% k2 = m2 - m1^2;
% k3 = m3 - 3*m2*m1 + 2*m1^3;
% k3_tild = k3 / (k2)^(3/2)
% % third feature
% m1 = mean(X(3, :));
% m2 = mean(X(3, :) .^ 2);
% m3 = mean(X(3, :) .^ 3);
% k2 = m2 - m1^2;
% k3 = m3 - 3*m2*m1 + 2*m1^3;
% k3_tild = k3 / (k2)^(3/2)
%% Covarience Matrix
mu_x = mean(X, 2);
Cx = (X - mu_x) * (X - mu_x)' / 1000;
[V, E] = eig(Cx);
[d, ind] = sort(diag(E), 'descend');
A = V(:, ind)';
y = A*X;
mu_y = mean(y, 2);
Cy = (y - mu_y) * (y - mu_y)' / 1000;
%% plot reduced diemension data in sensor space
Sel = 1;
X_red = A'*[y(1, :); zeros(2, 1000)];
figure;
scatter3(X_red(1, :), X_red(2, :), X_red(3, :))
xlim([-40 40])
ylim([-40 40])
zlim([-40 40])
xlabel('x1')
ylabel('x2')
zlabel('x3')
title('data after pca')
%% Visualize Data with Directions
scatter3(X(1, :), X(2, :), X(3, :))
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
hold on
p1 = mu_x';
p2 = p1+A(1, :)*sqrt(d(1));
plot3([p1(1), p2(1)],[p1(2), p2(2)],[p1(3), p2(3)], 'lineWidth', 2)
p2 = p1+A(2, :)*sqrt(d(2));
plot3([p1(1), p2(1)],[p1(2), p2(2)],[p1(3), p2(3)], 'lineWidth', 2)
p2 = p1+A(3, :)*sqrt(d(3));
plot3([p1(1), p2(1)],[p1(2), p2(2)],[p1(3), p2(3)], 'lineWidth', 2)
%% Covarience Matrix
mu_x = mean(X, 2);
Cx = (X - mu_x) * (X - mu_x)' / 1000;
[V, E] = eig(Cx);
[d, ind] = sort(diag(E), 'descend');
A = V(:, ind)';
B = sqrt(diag(d))^(-1);
D = B*A;
y = D*X;
mu_y = mean(y, 2);
Cy = (y - mu_y) * (y - mu_y)' / 1000;
scatter3(y(1, :), y(2, :), y(3, :))
xlabel('x1')
ylabel('x2')
zlabel('x3')
title('whitened and normalized data')
%% MATLAB PCA
[coeff, y1, latent] = pca(X');
mu_y = mean(y1, 1);
Cy = (y1 - mu_y)' * (y1 - mu_y) / 1000;
coeff
latent
%% SVD
N = size(X, 2);
[U, S, V] = svd(X);
U
D = S*S'/N

























