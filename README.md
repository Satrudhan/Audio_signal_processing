# Real-Time Sound Monitoring System

## Scenario
 Task is to develop a **real-time sound monitoring system** that detects abnormal frequency by analyzing audio signals from a microphone. The system should process the signal onboard, detect anomalies in frequency, and transmit data wirelessly to a cloud API for remote monitoring.

## Task Breakdown

### Part 1: Sensor Data Acquisition & Preprocessing
1. **Set up a microphone sensor** (Electret or MEMS) with the **ADC2** pin of an ESP32.
2. Sample the audio signal at **8 kHz or higher**. - Use the below mocked signal[0] or you can generate as well.
4. Mock Signal File: [0] [Link](https://github.com/Satrudhan/Audio_signal_processing/blob/main/InHive_Sensor/scripts/mock_signal.npy)
5. Implement **preprocessing techniques**:
   - Apply a **low-pass filter** (e.g., Butterworth) to remove high-frequency noise.
   - Normalize the signal.
   - Compute a **moving average** to smoothen variations.

### Part 2: Signal Processing & Anomaly Detection
1. Implement **Fast Fourier Transform (FFT)** to analyze frequency components.
2. Detect anomalies:
   - Identify **unexpected spikes or shifts** in frequency from 400Hz.
   - Use a **threshold-based approach** (e.g., if high-frequency noise exceeds X dB, trigger an alert).
   - Bonus: Use **Mel-Frequency Cepstral Coefficients (MFCCs)** for sound classification.

### Part 3: Wireless Data Transmission via Bluetooth & Wi-Fi
1. **InHive Sensor Firmware**:
   - Transmit **processed sound features using ESP WROOM 32** (not raw audio) over Bluetooth to the gateway.
   - Implement **data compression** for efficient transmission.
2. **Gateway Firmware**:
   - Receives data over BLE using ESP WROOM 32 from InHive Sensor.
   - Sends data to an **HTTP API endpoint**.
   - Implement **error handling** (e.g., retry failed transmissions).

## Implemented Parts
- **Firmware Code** (ESP32, C/C++ with Platform IO or ESP-IDF).
- **Python/Matlab scripts** for signal processing tests and analysis.

## Hardware Demo
[Hardware Demo Video](https://drive.google.com/file/d/1WkpvI5xJFhGXISHnsoeXTCJkSDkqH9Lf/view?usp=sharing)
