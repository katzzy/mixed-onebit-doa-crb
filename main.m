clear;
clc;
% Parameter settings

M1 = 5; % Number of subarrays 1
M2 = 5; % Number of subarrays 2
D = 2 * M2 * (M1 + 1) - 1; % Total number of virtual array sensors

theta = 30;
k = length(theta);
theta = sort(theta, 'ascend');
N = 100; % Number of snapshots
snr = -5;
non_quan_bits_list = [1, 2, 3, 4, 5, 6, 7, 8, 9]; % List of non-quantized array element counts
R = 1000; % Number of Monte Carlo tests

prompt = "Type:";
type = input(prompt);

switch type
    case 1
        disp("music SNR");
        SNR = [-10 -5  0 5 10 15 20];
        exp_SNR_analyze(M1, M2, k, N, SNR, theta, D, non_quan_bits_list, R);
    case 2
        disp("music N")
        N_values = [50 100 150 200 250 300 500];
        exp_N_analyze(M1, M2, k, N_values, snr, theta, D, non_quan_bits_list, R);
    case 3
        disp("music sensors")
        M1_values = [2,3,4,5,6];
        M2_values = [2,3,4,5,6];
        D_values = [11,23,39,59,83];
        exp_sensors_analyze(M1_values, M2_values, k, N, snr, theta, D_values, non_quan_bits_list, R);
    case 4
        disp("music angle")
        angle_set_values = -60:20:60;
        exp_angle_analyze(M1, M2, k, N, snr, angle_set_values, D, non_quan_bits_list, R);
    case 5
        disp("music undetermined")
        snr = 10;
        theta = [-54, -45, -37, -30, -22, -17, -11, -5, 4, 7, 15, 27, 31, 36, 45, 55];
        exp_undetermined_analyze(M1, M2, k, N, snr, theta, D, non_quan_bits_list, R);
    otherwise
        error('error type！');
end


% RMSE vs SNR 
function [] = exp_SNR_analyze(M1, M2, k, N, SNR, angle_set, D, non_quan_bits_list, R)

disp(angle_set)

rmae_results = zeros(length(non_quan_bits_list), length(SNR));
rmae_all_quan = zeros(1, length(SNR));
rmae_no_quan = zeros(1, length(SNR));

% Calculate the all-1-bit quantization benchmark
for j = 1:R
    for i = 1:length(SNR)
        angle_get_all_quan = music_experiment(M1, M2, k, N, SNR(i), angle_set, D, "all quan", 1);
        rmae_all_quan(i) = rmae_all_quan(i) + 1/k * sqrt(sum(abs(angle_set - angle_get_all_quan)));
    end
end
rmae_all_quan = rmae_all_quan/R;

% Calculate full-precision quantization benchmark
for j = 1:R
    for i = 1:length(SNR)
        angle_get_no_quan = music_experiment(M1, M2, k, N, SNR(i), angle_set, D, "no quan", 1);
        rmae_no_quan(i) = rmae_no_quan(i) + 1/k * sqrt(sum(abs(angle_set - angle_get_no_quan)));
    end
end
rmae_no_quan = rmae_no_quan/R;

% Calculate the results for different numbers of non-quantized array sensors
for nq_idx = 1:length(non_quan_bits_list)
    non_quan_bits = 1:non_quan_bits_list(nq_idx);
    rmae_temp = zeros(1, length(SNR));
    
    for j = 1:R
        for i = 1:length(SNR)
            angle_get_mix_quan = music_experiment(M1, M2, k, N, SNR(i), angle_set, D, "mix quan", non_quan_bits);
            rmae_temp(i) = rmae_temp(i) + 1/k * sqrt(sum(abs(angle_set - angle_get_mix_quan)));
        end
    end
    rmae_results(nq_idx, :) = rmae_temp/R;
end

% Drawing
colors = ['r', 'b', 'g', 'm', 'c','r', 'b', 'g', 'm'];
line_styles = {'-', '--', '-.', ':', '-', '--', '-.', ':', '-'};
figure(1);
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;

% Draw the all-1-bit quantization baseline
plot(SNR, rmae_all_quan, 'k--', 'LineWidth', 2, 'DisplayName', '1 bit for all');


% Plot the results for different numbers of unquantized array sensors
for nq_idx = 1:length(non_quan_bits_list)
    non_quan_bits = non_quan_bits_list(nq_idx);
    legend_str = sprintf('1 bit for %d sensors', M1 + M2 - non_quan_bits);
    
    plot(SNR, rmae_results(nq_idx, :), ...
         'Color', colors(nq_idx), ...
         'LineStyle', line_styles{nq_idx}, ...
         'LineWidth', 1.5, ...
         'DisplayName', legend_str);
end

% Draw the full-precision quantization baseline
plot(SNR, rmae_no_quan, 'y-.', 'LineWidth', 2, 'DisplayName', 'full-precision for all');

xlabel('SNR(dB)');
ylabel('RMSE (^\circ)');
legend('Location', 'best');
grid on;
title(num2str([M1 + M2, N, k], 'Array sensors: %d, snapshots: %d, sources: %d'));
hold off;

f = gcf;
exportgraphics(f, "SNR_comparison.pdf", "ContentType","vector")

end


% RMSE vs N
function [] = exp_N_analyze(M1, M2, k, N_values, snr, angle_set, D, non_quan_bits_list, R)



rmae_results = zeros(length(non_quan_bits_list), length(N_values));
rmae_all_quan = zeros(1, length(N_values));
rmae_no_quan = zeros(1, length(N_values));

% Calculate the all-1-bit quantization benchmark
for i = 1:length(N_values)
    for j = 1:R
        angle_get_all_quan = music_experiment(M1, M2, k, N_values(i), snr, angle_set, D, "all quan", non_quan_bits_list);
        rmae_all_quan(1, i) = rmae_all_quan(1, i) + 1/k * sqrt(sum(abs(angle_set - angle_get_all_quan)));
    end
    rmae_all_quan(1, i) = rmae_all_quan(1, i)/R;
end

% Calculate full-precision quantization benchmark
for i = 1:length(N_values)
    for j = 1:R
        angle_get_no_quan = music_experiment(M1, M2, k, N_values(i), snr, angle_set, D, "no quan", non_quan_bits_list);
        rmae_no_quan(1, i) = rmae_no_quan(1, i) + 1/k * sqrt(sum(abs(angle_set - angle_get_no_quan)));
    end
    rmae_no_quan(1, i) = rmae_no_quan(1, i)/R;
end



% Calculate the results for different numbers of non-quantized array sensors
for nq_idx = 1:length(non_quan_bits_list)
    non_quan_bits = 1:non_quan_bits_list(nq_idx);
    
    for i = 1:length(N_values)
        for j = 1:R
            angle_get_mix_quan = music_experiment(M1, M2, k, N_values(i), snr, angle_set, D, "mix quan", non_quan_bits);
            rmae_results(nq_idx, i) = rmae_results(nq_idx, i) + 1/k * sqrt(sum(abs(angle_set - angle_get_mix_quan)));
        end
        rmae_results(nq_idx, i) = rmae_results(nq_idx, i)/R;
    end
end

% 绘图
colors = ['r', 'b', 'g', 'm', 'c','r', 'b', 'g', 'm'];
line_styles = {'-', '--', '-.', ':', '-', '--', '-.', ':', '-'};
figure(1);
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;

% Draw the all-1-bit quantization baseline
plot(N_values, rmae_all_quan, 'k--', 'LineWidth', 2, 'DisplayName', '1 bit for all');

% Plot the results for different numbers of unquantized array sensors
for nq_idx = 1:length(non_quan_bits_list)
    non_quan_bits = non_quan_bits_list(nq_idx);
    legend_str = sprintf('1 bit for %d sensors', M1 + M2 - non_quan_bits);
    
    plot(N_values, rmae_results(nq_idx, :), ...
         'Color', colors(nq_idx), ...
         'LineStyle', line_styles{nq_idx}, ...
         'LineWidth', 1.5, ...
         'DisplayName', legend_str);
end

% Draw the full-precision quantization baseline
plot(N_values, rmae_no_quan, 'y-.', 'LineWidth', 2, 'DisplayName', 'full-precision for all');

xlabel('Snapshots');
ylabel('RMSE (^\circ)');
legend('Location', 'best');
grid on;
title(num2str([M1 + M2, snr, k], 'Array sensors: %d, SNR: %d, sources: %d'));
hold off;

f = gcf;
exportgraphics(f, "snapshot_comparison.pdf", "ContentType","vector")
end


% RMSE vs Sensors
function [] = exp_sensors_analyze(M1_values, M2_values, k, N, snr, angle_set, D_values, ~, R)



rmae_results = zeros(1, length(D_values));
rmae_all_quan = zeros(1, length(D_values));
rmae_no_quan = zeros(1, length(D_values));

% Calculate the all-1-bit quantization benchmark
for i = 1:length(D_values)
    for j = 1:R
        angle_get_all_quan = music_experiment(M1_values(i), M2_values(i), k, N, snr, angle_set, D_values(i), "all quan", 1);
        rmae_all_quan(1, i) = rmae_all_quan(1, i) + 1/k * sqrt(sum(abs(angle_set - angle_get_all_quan)));
    end
    rmae_all_quan(1, i) = rmae_all_quan(1, i)/R;
end

% Calculate full-precision quantization benchmark
for i = 1:length(D_values)
    for j = 1:R
        angle_get_no_quan = music_experiment(M1_values(i), M2_values(i), k, N, snr, angle_set, D_values(i), "no quan", 1);
        rmae_no_quan(1, i) = rmae_no_quan(1, i) + 1/k * sqrt(sum(abs(angle_set - angle_get_no_quan)));
    end
    rmae_no_quan(1, i) = rmae_no_quan(1, i)/R;
end


non_quan_bits = 1;

for i = 1:length(D_values)
    
    for j = 1:R
        angle_get_mix_quan = music_experiment(M1_values(i), M2_values(i), k, N, snr, angle_set, D_values(i), "mix quan", non_quan_bits);
        rmae_results(1, i) = rmae_results(1, i) + 1/k * sqrt(sum(abs(angle_set - angle_get_mix_quan)));
    end
    rmae_results(1, i) = rmae_results(1, i)/R;
end


% Drawing
figure(1);
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;

plot(D_values, rmae_all_quan, 'k--', 'LineWidth', 2, 'DisplayName', '1 bit for all');

plot(D_values, rmae_results, 'r-', 'LineWidth', 2, 'DisplayName', '1 bit for 1 sensors');

plot(D_values, rmae_no_quan, 'y-.', 'LineWidth', 2, 'DisplayName', 'full-precision for all');


xlabel('Sensors');
ylabel('RMSE (^\circ)');
legend('Location', 'best');
grid on;
title(num2str([snr, N, k], 'SNR: %d, snapshots: %d, sources: %d'));
hold off;

f = gcf;
exportgraphics(f, "elements_comparison.pdf", "ContentType","vector")
end

% RMSE vs Angle
function [] = exp_angle_analyze(M1, M2, k, N, snr, angle_set_values, D, non_quan_bits_list, R)



rmae_results = zeros(length(non_quan_bits_list), length(angle_set_values));
rmae_all_quan = zeros(1, length(angle_set_values));
rmae_no_quan = zeros(1, length(angle_set_values));

% Calculate the all-1-bit quantization benchmark
for i = 1:length(angle_set_values)
    for j = 1:R
        angle_get_all_quan = music_experiment(M1, M2, k, N, snr, angle_set_values(i), D, "all quan", 1);
        rmae_all_quan(1, i) = rmae_all_quan(1, i) + 1/k * sqrt(sum(abs(angle_set_values(i) - angle_get_all_quan)));
    end
    rmae_all_quan(1, i) = rmae_all_quan(1, i)/R;
end

for i = 1:length(angle_set_values)
    for j = 1:R
        angle_get_no_quan = music_experiment(M1, M2, k, N, snr, angle_set_values(i), D, "no quan", 1);
        rmae_no_quan(1, i) = rmae_no_quan(1, i) + 1/k * sqrt(sum(abs(angle_set_values(i) - angle_get_no_quan)));
    end
    rmae_no_quan(1, i) = rmae_no_quan(1, i)/R;
end

% Calculate the results for different numbers of non-quantized array sensors
for nq_idx = 1:length(non_quan_bits_list)
    non_quan_bits = 1:non_quan_bits_list(nq_idx);
    
    for i = 1:length(angle_set_values)
        for j = 1:R
            angle_get_mix_quan = music_experiment(M1, M2, k, N, snr, angle_set_values(i), D, "mix quan", non_quan_bits);
            rmae_results(nq_idx, i) = rmae_results(nq_idx, i) + 1/k * sqrt(sum(abs(angle_set_values(i) - angle_get_mix_quan)));
        end
        rmae_results(nq_idx, i) = rmae_results(nq_idx, i)/R;
    end
end

% Drawing
colors = ['r', 'b', 'g', 'm', 'c','r', 'b', 'g', 'm'];
line_styles = {'-', '--', '-.', ':', '-', '--', '-.', ':', '-'};
figure(1);
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;

% Draw the full-precision quantization baseline
plot(angle_set_values, rmae_all_quan, 'k--', 'LineWidth', 2, 'DisplayName', '1 bit for all');

% Plot the results for different numbers of unquantized array sensors
for nq_idx = 1:length(non_quan_bits_list)
    non_quan_bits = non_quan_bits_list(nq_idx);
    legend_str = sprintf('1 bit for %d sensors', M1 + M2 - non_quan_bits);

    
    plot(angle_set_values, rmae_results(nq_idx, :), ...
         'Color', colors(nq_idx), ...
         'LineStyle', line_styles{nq_idx}, ...
         'LineWidth', 1.5, ...
         'DisplayName', legend_str);
end

plot(angle_set_values, rmae_no_quan, 'y-.', 'LineWidth', 2, 'DisplayName', 'full-precision for all');

xlabel('Incident Signal Angle (^\circ)');
ylabel('RMSE (^\circ)');
legend('Location', 'north');
grid on;
title(num2str([M1 + M2, snr,N, k], 'Array sensors: %d, SNR: %d, snapshots: %d, sources: %d'));
hold off;

f = gcf;
exportgraphics(f, "angle_comparison.pdf", "ContentType","vector")
end


% undetermind
function [] = exp_undetermined_analyze(M1, M2, ~, N, snr, theta, D, non_quan_bits_list, R)
angle_set_number = [10, 11, 12, 13, 14, 15, 16];


rmae_results = zeros(length(non_quan_bits_list), length(angle_set_number));
rmae_all_quan = zeros(1, length(angle_set_number));
rmae_no_quan = zeros(1, length(angle_set_number));

% Calculate the all-1-bit quantization benchmark
for j = 1:R
    for i = angle_set_number
        angle_set = theta(1:i);
        k = length(angle_set);
        angle_set = sort(angle_set, 'ascend');
        angle_get_all_quan = music_experiment(M1, M2, k, N, snr, angle_set, D, "all quan", 1);
        rmae_all_quan(i - 9) = rmae_all_quan(i - 9) + 1/k * sqrt(sum(abs(angle_set - angle_get_all_quan)));
    end
end
rmae_all_quan = rmae_all_quan/R;

% Calculate full-precision quantization benchmark
for j = 1:R
    for i = angle_set_number
        angle_set = theta(1:i);
        k = length(angle_set);
        angle_set = sort(angle_set, 'ascend');
        angle_get_no_quan = music_experiment(M1, M2, k, N, snr, angle_set, D, "no quan", 1);
        rmae_no_quan(i - 9) = rmae_no_quan(i - 9) + 1/k * sqrt(sum(abs(angle_set - angle_get_no_quan)));
    end
end
rmae_no_quan = rmae_no_quan/R;

% Calculate the results for different numbers of non-quantized array sensors
for nq_idx = 1:length(non_quan_bits_list)
    non_quan_bits = 1:non_quan_bits_list(nq_idx);
    rmae_temp = zeros(1, length(angle_set_number));

    for j = 1:R
        for i = angle_set_number
            angle_set = theta(1:i);
            k = length(angle_set);
            angle_set = sort(angle_set, 'ascend');
            angle_get_mix_quan = music_experiment(M1, M2, k, N, snr, angle_set, D, "mix quan", non_quan_bits);
            rmae_temp(i - 9) = rmae_temp(i - 9) + 1/k * sqrt(sum(abs(angle_set - angle_get_mix_quan)));
        end
    end
    rmae_results(nq_idx, :) = rmae_temp/R;
end

% Drawing
colors = ['r', 'b', 'g', 'm', 'c','r', 'b', 'g', 'm'];
line_styles = {'-', '--', '-.', ':', '-', '--', '-.', ':', '-'};
figure(1);
set(gcf, 'Position', [100, 100, 1000, 900]);
hold on;

% Draw the all-1-bit quantization baseline
plot(angle_set_number, rmae_all_quan, 'k--', 'LineWidth', 2, 'DisplayName', '1 bit for all');


% Plot the results for different numbers of unquantized array sensors
for nq_idx = 1:length(non_quan_bits_list)
    non_quan_bits = non_quan_bits_list(nq_idx);
    legend_str = sprintf('1 bit for %d sensors', M1 + M2 - non_quan_bits);

    plot(angle_set_number, rmae_results(nq_idx, :), ...
         'Color', colors(nq_idx), ...
         'LineStyle', line_styles{nq_idx}, ...
         'LineWidth', 1.5, ...
         'DisplayName', legend_str);
end

% Draw the full-precision quantization baseline
plot(angle_set_number, rmae_no_quan, 'y-.', 'LineWidth', 2, 'DisplayName', 'full-precision for all');

xlabel('Angle Number');
ylabel('RMSE (^\circ)');
legend('Location', 'northwest');
grid on;
title(num2str([M1 + M2, snr, N], 'Array sensors: %d, SNR: %d, snapshots: %d'));
hold off;

f = gcf;
exportgraphics(f, "undetermined.pdf", "ContentType","vector")

end