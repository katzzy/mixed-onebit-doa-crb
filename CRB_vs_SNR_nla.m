N = 100;        %Number of snapshots

% M = 15;         %Number of sensors
% D = 0:M-1;      % Position efficients of sensors
% 
K = 3;          %Number of sources
theta = [-21 2 19];     % Impinging angles of sources
omega = [-1i -3i 2i];   % Exponents of complex exponential source signals

M1 = 5;         % -- NLA --
M2 = 5;         % -- sensors --
M = M1 + M2;
D = [0:M1-1,(1:M2)*(M1+1)-1];

A = exp(-1i*pi*D'*sind(theta));  % steering matrix
dA = -1i*pi*D'*cosd(theta).*A;   % derevatives of steering matrix
S = exp((0:N-1)'*omega);
epsilon = S(:);                 % vectorized vec(S)
B = kron(eye(N), A);            % The observation matrix after vectorizing S

Delta = zeros(N*M, K);
for k = 1:K
    Delta(:,k) = kron(S(:,k), dA(:,k));
end

Lambda = [real(Delta); imag(Delta)];
C = [real(B), -imag(B); imag(B), real(B)];
u = [real(epsilon); imag(epsilon)];

SNR = -10:2:20; % Signal to Noise Ratio in dB
crbx = zeros(M+1, size(SNR,2));
H = sqrt(2) * C * u; % rearranged signal samples
SE = H' * H / ( 2 * M * N ); % average signal power
for k = 1:size(SNR,2)
    sigma = sqrt( SE * 10^(-SNR(k) / 10) );
    h = H / sigma;
    d1 = normcdf(h).*normcdf(-h)./(normpdf(h).^2);
    m = 0;  % number of sensors left for full precision quantization
%% CRB  
    while m <= M
        E = diag(sqrt(d1.^-1)) * Lambda;
        F = diag(sqrt(d1.^-1)) * C;
        Fc = eye(2*M*N) - F*pinv(F);    % the orthognal complement matrix of F
        CRB = sigma^2 / 2 * inv(E' * Fc * E);
        crb = min(sqrt(diag(CRB))) * 180 / pi;
        m = m + 1;
        crbx(m, k) = crb;
        d1(m:M:end) = 1;    % set the values corresponding to 1:m sensors to be 1, on those sensors 1bit quantization is not applied
    end 
        
    display([SNR(k) crbx(:,k)'])
end

figure('Position', [100, 100, 1200, 800]);  % Adjust figure size here
hold on;
semilogy(SNR, crbx(1,:), '-s', repmat(SNR', 1, M-1), crbx(2:M,:)', SNR, crbx(M+1,:), '-d'), grid on
xlim([SNR(1), SNR(end)]);
set(gca, 'XTick', SNR, 'FontSize', 24);
legend(['1bit for all'; cellstr(num2str([M-1:-1:1]', '1bit for %02d sensors')); ' full-precision quantization'], 'FontSize', 22)
ylabel('$\sqrt {CRB} (^\circ)$', 'Interpreter', 'latex', 'FontSize', 28)
xlabel('SNR(dB)', 'FontSize', 28)
title(num2str([M, N, K], 'Array sensors: %d, snapshots: %d, sources: %d'), 'FontSize', 30);
hold off;
f = gcf;
exportgraphics(f, "CRB_NLA.pdf", "ContentType","vector")
