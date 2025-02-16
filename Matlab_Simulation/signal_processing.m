% Real-Time Signal Processing with FFT & Filtering
% Author: Satrudhan Kumar Chauhan
% Description: This script processes a mock signal in real-time by normalizing,
% filtering, computing FFT, and plotting results dynamically.

clc; clear; close all;

% Load mock signal from a NumPy file
mock_signal = readNPY('mock_signal.npy'); 

% =================== Signal Processing Parameters ===================
SIGNAL_SIZE = 256;  % Number of samples per processing chunk
SAMPLE_RATE = 8000; % Sampling frequency in Hz
CUTOFF_FREQ = 1000; % Low-pass filter cutoff frequency (1kHz)
FILTER_ORDER = 3;   % Butterworth filter order

% Ensure mock_signal is a row vector for processing
if size(mock_signal, 1) > 1
    mock_signal = mock_signal';
end

% =================== Design Low-Pass Butterworth Filter ===================
% A Butterworth filter smooths high-frequency noise from the signal
[b, a] = butter(FILTER_ORDER, CUTOFF_FREQ / (SAMPLE_RATE / 2), 'low');

% =================== Initialize Processing Variables ===================
mockSignalIndex = 1;              % Start index for reading the signal
totalSamples = length(mock_signal); % Total number of samples in the signal

% =================== Create Real-Time Plot ===================
figure;

% Subplot 1: Original vs Normalized Signal
subplot(2,2,2);
originalVsNormalizedPlot = plot(zeros(1, SIGNAL_SIZE), 'k', 'DisplayName', 'Original');
hold on;
normalizedPlot = plot(zeros(1, SIGNAL_SIZE), 'b', 'DisplayName', 'Normalized');
hold off;
xlabel('Samples'); ylabel('Amplitude'); title('Original vs Normalized Signal');
legend;

% Subplot 2: Original vs Filtered Signal
subplot(2,2,1);
originalVsFilteredPlot = plot(zeros(1, SIGNAL_SIZE), 'k', 'DisplayName', 'Original');
hold on;
filteredPlot = plot(zeros(1, SIGNAL_SIZE), 'g', 'DisplayName', 'Filtered');
hold off;
xlabel('Samples'); ylabel('Amplitude'); title('Original vs Filtered Signal');
legend;

% Subplot 3: Frequency Spectrum (FFT) - Now includes all signals
subplot(2,2,4);
fftPlot_original = plot(zeros(1, SIGNAL_SIZE/2+1), 'k', 'DisplayName', 'Original Spectrum');
hold on;
fftPlot_filtered = plot(zeros(1, SIGNAL_SIZE/2+1), 'g', 'DisplayName', 'Filtered Spectrum');
fftPlot_smoothed = plot(zeros(1, SIGNAL_SIZE/2+1), 'm', 'DisplayName', 'Smoothed Spectrum');
hold off;
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); title('FFT Spectrum');
xlim([0 SAMPLE_RATE/2]); 
legend;

% Subplot 4: Smoothed Signal
subplot(2,2,3);
smoothPlot = plot(zeros(1, SIGNAL_SIZE), 'm');
xlabel('Samples'); ylabel('Amplitude'); title('Smoothed Signal');

%% =================== Real-Time Signal Processing Loop ===================
while mockSignalIndex + SIGNAL_SIZE - 1 <= totalSamples
    % 1. Read a chunk of the signal (256 samples at a time)
    original_signal = mock_signal(mockSignalIndex:mockSignalIndex+SIGNAL_SIZE-1);
    mockSignalIndex = mockSignalIndex + SIGNAL_SIZE; % Move to the next chunk

    % 2. Normalize the signal to a range of [-1, 1]
    normalized_signal = normalizeSignal(original_signal);

    % 3. Apply a low-pass Butterworth filter to remove noise
    filtered_signal = filtfilt(b, a, normalized_signal); 

    % 4. Smooth the signal using a moving average filter
    smoothed_signal = smoothdata(filtered_signal, 'movmean', 10);

    % 5. Compute FFT for frequency domain analysis
    [freqBins, fftOriginal] = computeFFT(original_signal, SAMPLE_RATE);
    [~, fftFiltered] = computeFFT(filtered_signal, SAMPLE_RATE);
    [~, fftSmoothed] = computeFFT(smoothed_signal, SAMPLE_RATE);

    % =================== Update Plots in Real-Time ===================
    set(originalVsNormalizedPlot, 'YData', original_signal);
    set(normalizedPlot, 'YData', normalized_signal);
    
    set(originalVsFilteredPlot, 'YData', original_signal);
    set(filteredPlot, 'YData', filtered_signal);
    
    % Update FFT plots for all three signals
    set(fftPlot_original, 'XData', freqBins, 'YData', fftOriginal);
    set(fftPlot_filtered, 'XData', freqBins, 'YData', fftFiltered);
    set(fftPlot_smoothed, 'XData', freqBins, 'YData', fftSmoothed);

    set(smoothPlot, 'YData', smoothed_signal);

    pause(0.1); % Simulates real-time processing delay
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
    N = length(signal); % Number of samples
    Y = fft(signal);    % Compute Fast Fourier Transform
    P2 = abs(Y / N);    % Normalize the FFT result
    P1 = P2(1:N/2+1);   % Take the positive frequency spectrum
    P1(2:end-1) = 2 * P1(2:end-1); % Scale appropriately
    fftMagnitudes = 20 * log10(P1 + eps); % Convert magnitude to dB scale
    freqBins = sample_rate * (0:(N/2)) / N; % Generate corresponding frequency bins
end