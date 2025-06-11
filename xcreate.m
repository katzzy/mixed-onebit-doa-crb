function [x_normalization]=xcreate(M1,M2,k,N,snr,angle_set)

d = zeros(1, M1 + M2);
for i = 1:M1+M2
    if i <= M1
        d(i) = i - 1;
    else
        d(i) = M1 + (i - M1 - 1)*(M1 + 1);
    end
end

c = zeros(1,k); % Energy of each source signal s(t)
while min(c)<0.4 % Ensure that the energy of each randomly generated source signal s(t) is greater than 0.4
    s = randn(k, N);
    % Calculate energy
    for i = 1:k
        c(i) = sum(s(i,:).^2)/length(s(i,:));
    end
end

% R_corr = corrcoef(s'); 
% disp('Correlation coefficient matrix:');
% disp(abs(R_corr));

% Generate the direction vector a(Θ)
a = zeros(M1 + M2, k);
for j = 1:k
    a(:, j) = exp(-1i * pi * sin(angle_set(j)) * d');
end

% Generate noise based on signal-to-noise ratio and produce output
x = a*s;
% Load Gaussian white noise
x_noisy = awgn(x, snr, 'measured');

% Zero-mean normalization
Ex = mean(x_noisy, 2);
Ex_repeat = repmat(Ex, 1, N);
x_noisy_zero_mean = x_noisy - Ex_repeat;

x_normalization = normalization_x(x_noisy_zero_mean, N);

end


function [x_normal] = normalization_x(x, N)
    D_real = mean(real(x).^2,2);
    D_imag = mean(imag(x).^2,2);
    Sx = sqrt(D_imag + D_real);
    Dx_repeat = repmat(Sx, 1, N);
    x_normal = x./Dx_repeat;
end