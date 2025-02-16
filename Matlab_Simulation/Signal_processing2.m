% Real-Time Signal Processing with FFT & Anomaly Detection
% Author: Satrudhan Kumar Chauhan
% Description: This script processes a mock signal in real-time by normalizing,
% filtering, computing FFT, detecting anomalies, and plotting results dynamically.

clc; clear; close all;

% =================== Load Mock Signal ===================
mock_signal = readNPY('mock_signal.npy'); 

% =================== Signal Processing Parameters ===================
SIGNAL_SIZE = 256;  % Number of samples per processing chunk
SAMPLE_RATE = 8000; % Sampling frequency in Hz
CUTOFF_FREQ = 1000; % Low-pass filter cutoff frequency (1kHz)
FILTER_ORDER = 3;   % Butterworth filter order
ANOMALY_FREQ = 400; % Frequency to monitor for anomalies
ANOMALY_THRESHOLD = -20; % Threshold in dB for anomaly detection

% Ensure mock_signal is a row vector for processing
if size(mock_signal, 1) > 1
    mock_signal = mock_signal';
end

% =================== Design Low-Pass Butterworth Filter ===================
[b, a] = butter(FILTER_ORDER, CUTOFF_FREQ / (SAMPLE_RATE / 2), 'low');

% =================== Initialize Processing Variables ===================
mockSignalIndex = 1;
totalSamples = length(mock_signal);

% =================== Open Figures for Real-Time Visualization ===================
figure('Name', 'Time-Domain Analysis');
subplot(3,1,1); title('Original Signal'); hold on; grid on;
subplot(3,1,2); title('Filtered Signal'); hold on; grid on;
subplot(3,1,3); title('Smoothed Signal'); hold on; grid on;

figure('Name', 'Frequency-Domain Analysis');
subplot(3,1,1); title('FFT - Original Signal'); hold on; grid on;
subplot(3,1,2); title('FFT - Filtered Signal'); hold on; grid on;
subplot(3,1,3); title('FFT - Smoothed Signal'); hold on; grid on;

% =================== Real-Time Signal Processing Loop ===================
while mockSignalIndex + SIGNAL_SIZE - 1 <= totalSamples
    % 1. Read a chunk of the signal
    original_signal = mock_signal(mockSignalIndex:mockSignalIndex+SIGNAL_SIZE-1);
    mockSignalIndex = mockSignalIndex + SIGNAL_SIZE;

    % 2. Normalize the signal
    normalized_signal = normalizeSignal(original_signal);

    % 3. Apply low-pass Butterworth filter
    filtered_signal = filtfilt(b, a, normalized_signal);

    % 4. Smooth the signal using a moving average filter
    smoothed_signal = smoothdata(filtered_signal, 'movmean', 10);

    % 5. Compute FFT
    [freqBins, fftOriginal] = computeFFT(original_signal, SAMPLE_RATE);
    [~, fftFiltered] = computeFFT(filtered_signal, SAMPLE_RATE);
    [~, fftSmoothed] = computeFFT(smoothed_signal, SAMPLE_RATE);

    % =================== Plot Time-Domain Signals ===================
    figure(1);
    subplot(3,1,1); cla; plot(original_signal, 'r'); title('Original Signal');
    xlabel('Samples'); ylabel('Amplitude'); grid on;

    subplot(3,1,2); cla; plot(filtered_signal, 'b'); title('Filtered Signal');
    xlabel('Samples'); ylabel('Amplitude'); grid on;

    subplot(3,1,3); cla; plot(smoothed_signal, 'g'); title('Smoothed Signal');
    xlabel('Samples'); ylabel('Amplitude'); grid on;

    % =================== Plot Frequency-Domain Signals ===================
    figure(2);
    subplot(3,1,1); cla; plot(freqBins, fftOriginal, 'r'); title('FFT - Original Signal');
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on;

    subplot(3,1,2); cla; plot(freqBins, fftFiltered, 'b'); title('FFT - Filtered Signal');
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on;

    subplot(3,1,3); cla; plot(freqBins, fftSmoothed, 'g'); title('FFT - Smoothed Signal');
    xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on;

    % =================== Anomaly Detection ===================
    [~, idx] = min(abs(freqBins - ANOMALY_FREQ));
    if fftOriginal(idx) > ANOMALY_THRESHOLD
        fprintf('Anomaly detected at 400Hz! Magnitude: %.2f dB\n', fftOriginal(idx));

        % Mark anomaly on FFT plot
        figure(2);
        subplot(3,1,1);
        hold on; plot(freqBins(idx), fftOriginal(idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    end
    % =================== Anomaly Detection ===================
    % [~, idx] = min(abs(freqBins - ANOMALY_FREQ));
    % if fftFiltered(idx) > ANOMALY_THRESHOLD
    %     fprintf('Anomaly detected at 400Hz! Magnitude: %.2f dB\n', fftFiltered(idx));
    % 
    %     % Mark anomaly on FFT - Filtered Signal plot
    %     figure(2);
    %     subplot(3,1,2);
    %     hold on; plot(freqBins(idx), fftFiltered(idx), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    % end
    % 
    
    drawnow; % Update plots in real-time
    pause(0.1); % Simulate real-time delay
end



%% =================== Helper Functions ===================

% Function to normalize signal to range [-1, 1]
function normalized_signal = normalizeSignal(signal)
    min_val = min(signal);
    max_val = max(signal);
    range = max_val - min_val;
    if range == 0
        normalized_signal = signal;
    else
        normalized_signal = 2 * (signal - min_val) / range - 1;
    end
end

% Function to compute FFT and return frequency bins and magnitudes
function [freqBins, fftMagnitudes] = computeFFT(signal, sample_rate)
    N = length(signal);
    Y = fft(signal);
    P2 = abs(Y / N);
    P1 = P2(1:N/2+1);
    P1(2:end-1) = 2 * P1(2:end-1);
    fftMagnitudes = 20 * log10(P1 + eps); % Convert to dB
    freqBins = sample_rate * (0:(N/2)) / N;
end
