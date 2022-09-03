function y = smoother(x, a_thr, t_thr, high_freq, feq)
y = x > a_thr;
tmp = zeros(size(y));
L = t_thr*feq;
for i = 1:length(x)
    if i - L / 2 <= 0
        first_half = y(1:i-1);
        second_half = y(i+1:i+L/2);
    elseif i + L / 2 > length(x)
        first_half = y(i-L/2:i-1);
        second_half = y(i+1:end);
    else
        first_half = y(i-L/2:i-1);
        second_half = y(i+1:i+L/2);
    end
    chunk = [first_half second_half];
    if sum(chunk) > high_freq
        tmp(i) = 1;
    else
        tmp(i) = 0;
    end
end
y = tmp;