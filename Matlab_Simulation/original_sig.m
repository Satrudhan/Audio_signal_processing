clc; clear; close all;
%% User-defined number of samples
num_samples = 256;  % Adjust this value as needed

%% Load the signal
signal_full = readNPY('mock_signal.npy');  % Load full signal
fs = 8000;  % Sampling rate (8 kHz)

% Ensure we do not exceed available data
num_samples = min(num_samples, length(signal_full));
signal = signal_full(1:num_samples);  

t = (0:length(signal)-1) / fs;  % Time vector

%% Design and Apply the Butterworth Low-Pass Filter
cutoff_freq = 1000;  % Cutoff frequency (1 kHz)
order = 3;  % 3rd-order Butterworth filter

% Butterworth filter design
[b, a] = butter(order, cutoff_freq/(fs/2), 'low');

% Apply the filter using zero-phase filtering
filtered_signal = filtfilt(b, a, signal);

%% Normalize Both Signals to [-1,1]
normalized_original = signal / max(abs(signal));  
normalized_filtered = filtered_signal / max(abs(filtered_signal));  

%% Apply Moving Average Filter to the Filtered Signal
window_size = 50;  % Number of samples for averaging
moving_avg_filter = ones(1, window_size) / window_size;  % Create filter

smoothed_signal = conv(filtered_signal, moving_avg_filter, 'same');  % Apply filter

%% Normalize the Smoothed Signal
normalized_smoothed = smoothed_signal / max(abs(smoothed_signal));

%% FFT Analysis for Frequency Components
N = length(signal);  % Number of samples
f = (0:N-1) * (fs/N);  % Frequency vector

% Compute FFT for original, filtered, and smoothed signals
fft_original = abs(fft(signal));
fft_filtered = abs(fft(filtered_signal));
fft_smoothed = abs(fft(smoothed_signal));

% Convert FFT to dB scale
fft_original_dB = 20*log10(fft_original);
fft_filtered_dB = 20*log10(fft_filtered);
fft_smoothed_dB = 20*log10(fft_smoothed);

%% Plot all graphs in a single figure
figure;

% 1. Original vs Filtered Signal
subplot(2,2,1);
plot(t, signal, 'r', 'DisplayName', 'Original Signal');
hold on;
plot(t, filtered_signal, 'b', 'DisplayName', 'Filtered Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Original vs Filtered Signal');
legend; grid on;
hold off;

% 2. Original vs Normalized Signal
subplot(2,2,2);
plot(t, normalized_original, 'r', 'DisplayName', 'Normalized Original Signal');
hold on;
plot(t, normalized_filtered, 'b', 'DisplayName', 'Normalized Filtered Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Original vs Normalized Signal');
legend; grid on;
hold off;

% 3. FFT Spectrum (Original, Filtered, Smoothed) in a single subplot
subplot(2,2,3);
plot(f(1:floor(N/2)), fft_original_dB(1:floor(N/2)), 'r', 'LineWidth', 1.5, 'DisplayName', 'Original FFT');
hold on;
plot(f(1:floor(N/2)), fft_filtered_dB(1:floor(N/2)), 'b', 'LineWidth', 1.5, 'DisplayName', 'Filtered FFT');
plot(f(1:floor(N/2)), fft_smoothed_dB(1:floor(N/2)), 'g', 'LineWidth', 1.5, 'DisplayName', 'Smoothed FFT');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('FFT Spectrum (Original, Filtered, Smoothed)');
legend; grid on;
hold off;

% 4. Smoothed Signal
subplot(2,2,4);
plot(t, smoothed_signal, 'g', 'DisplayName', 'Smoothed Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Smoothed Signal');
legend; grid on;

% Adjust figure size and spacing
sgtitle('Signal Processing Analysis');  % Set overall title
set(gcf, 'Position', [100, 100, 1000, 600]); % Adjust figure size

fprintf('All plots generated successfully in a single figure.\n');
