SNR = 5; % Fixed Signal to Noise Ratio in dB

N_values = 50:50:500; % Range of snapshots from 50 to 500

K = 3; % Number of sources
theta = [-21 2 19]; % Impinging angles of sources
omega = [-1i -3i 2i]; % Exponents of complex exponential source signals

M = 10;
D = 0:M-1;

crbx = zeros(M+1, length(N_values));

for n_idx = 1:length(N_values)
    N = N_values(n_idx);
    
    A = exp(-1i*pi*D'*sind(theta)); % steering matrix
    dA = -1i*pi*D'*cosd(theta).*A; % derivatives of steering matrix
    S = exp((0:N-1)'*omega);
    epsilon = S(:); % vectorized vec(S)
    B = kron(eye(N), A); % The observation matrix after vectorizing S

    Delta = zeros(N*M, K);
    for k = 1:K
        Delta(:,k) = kron(S(:,k), dA(:,k));
    end

    Lambda = [real(Delta); imag(Delta)];
    C = [real(B), -imag(B); imag(B), real(B)];
    u = [real(epsilon); imag(epsilon)];

    H = sqrt(2) * C * u; % rearranged signal samples
    SE = H' * H / ( 2 * M * N ); % average signal power

    sigma = sqrt( SE * 10^(-SNR / 10) );
    h = H / sigma;
    d1 = normcdf(h).*normcdf(-h)./(normpdf(h).^2);
    m = 0;  % number of sensors left for full precision quantization
    
    %% CRB  
    while m <= M
        E = diag(sqrt(d1.^-1)) * Lambda;
        F = diag(sqrt(d1.^-1)) * C;
        Fc = eye(2*M*N) - F*pinv(F);    % the orthogonal complement matrix of F
        CRB = sigma^2 / 2 * inv(E' * Fc * E);
        crb = min(sqrt(diag(CRB))) * 180 / pi;
        m = m + 1;
        crbx(m, n_idx) = crb;
        d1(m:M:end) = 1;    % set the values corresponding to 1:m sensors to be 1, on those sensors 1bit quantization is not applied
    end 
end

figure('Position', [100, 100, 1200, 800]);  % Adjust figure size here
hold on;
plot(N_values, crbx(1,:), '-s', 'LineWidth', 1.5);
for m = 2:M
    plot(N_values, crbx(m,:), 'LineWidth', 1.5);
end
plot(N_values, crbx(M+1,:), '-d', 'LineWidth', 1.5);
grid on;
xlim([N_values(1), N_values(end)]);
set(gca, 'XTick', N_values);
legend(['1bit for all'; cellstr(num2str([M-1:-1:1]', '1bit for %02d sensors')); ' full-precision quantization'], 'FontSize', 14)
ylabel('$\sqrt {CRB} (^\circ)$', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('Snapshots', 'FontSize', 16)
title(num2str([M, K, SNR], 'Array sensors: %d, sources: %d, SNR: %d dB'), 'FontSize', 18);
hold off;
f = gcf;
exportgraphics(f, "CRB_vs_N_ULA.pdf", "ContentType","vector")