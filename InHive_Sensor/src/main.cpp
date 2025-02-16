/* InHive_Sensor.ino */
#include <Arduino.h>
#include <arduinoFFT.h>
#include <BLEDevice.h>
#include <BLEServer.h>
#include <BLEUtils.h>
#include <BLE2902.h>
#include "mock_signal.h" // Include mock signal header

// ------------------ Signal Processing Constants ------------------
#define SIGNAL_SIZE         256
#define SAMPLE_RATE         8000
#define TARGET_FREQUENCY    400
#define THRESHOLD_DB        10
#define UPDATE_INTERVAL     100

// ------------------ Additional Feature Extraction ---------------
#define TOP_PEAKS_COUNT     3

// ------------------ BLE Constants -------------------------------
#define SERVICE_UUID        "91bad492-b950-4226-aa2b-4ede9fa42f59"
#define CHARACTERISTIC_UUID "cba1d466-344c-4be3-ab3f-189f80dd7518"

// ------------------ Data structure for sound features -----------
struct SoundFeaturesCompressed {
    uint16_t dominantFrequency;  // 2 bytes
    uint16_t magnitude;          // 2 bytes
    uint16_t frequencyDeviation; // 2 bytes
    // uint32_t timestamp;          // 4 bytes
    uint16_t topFrequency1;      // 2 bytes
    uint16_t topFrequency2;      // 2 bytes
    uint16_t topFrequency3;      // 2 bytes
    uint8_t topMagnitude1;       // 1 byte
    uint8_t topMagnitude2;       // 1 byte
    uint8_t topMagnitude3;       // 1 byte
    uint8_t anomalyCount;        // 1 byte
};  // *Total: 20 bytes*
struct SoundFeatures {
    float  dominantFrequency;
    float  magnitude;
    float  frequencyDeviation;
    // uint32_t timestamp;          // keep 4-byte fields together
    float  topFrequency1;
    float  topFrequency2;
    float  topFrequency3;
    float  topMagnitude1;
    float  topMagnitude2;
    float  topMagnitude3;
    float  psd;
    float  rmsValue;
    float  zcr;
    uint8_t anomalyCount;        // place the 1-byte field at the end
};

// ------------------ Global Variables -----------------------------
ArduinoFFT<float> FFT;
BLEServer* pServer = nullptr;
BLECharacteristic* pCharacteristic = nullptr;
bool deviceConnected = false;

float* processed_signal;   
float* vReal;             
float* vImag;             
unsigned long lastUpdateTime = 0;
size_t mockSignalIndex = 0;  // Index to track position in mock signal

// Simple circular buffer for signal processing
template <typename T, size_t S>
class SimpleBuffer {
    T buffer[S];
    size_t head = 0;
    size_t count = 0;

public:
    void push(T value) {
        buffer[head] = value;
        head = (head + 1) % S;
        if (count < S) count++;
    }

    T get(size_t index) const {
        return buffer[(S + head - count + index) % S];
    }

    size_t size() const { return count; }
    void clear() { head = 0; count = 0; }
};

// Statistics tracking
struct AnomalyStats {
    int anomalyCount = 0;
    float avgFrequencyDeviation = 0;
    float maxMagnitude = 0;
    unsigned long lastAnomalyTime = 0;
};

AnomalyStats stats;
SimpleBuffer<float, 16> signalBuffer;

// ------------------ BLE Server Callbacks -------------------------
class ServerCallbacks: public BLEServerCallbacks {
    void onConnect(BLEServer* pServer) override {
        deviceConnected = true;
        Serial.println("Gateway connected!");
    }
    
    void onDisconnect(BLEServer* pServer) override {
        deviceConnected = false;
        Serial.println("Gateway disconnected!");
        pServer->getAdvertising()->start();
    }
};

// ------------------ Utility Functions ----------------------------
// Return +1 if x>0, -1 if x<0, 0 if x=0
int signOf(float x) {
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

// ------------------ Signal Processing Functions ------------------
// Normalize signal to range [-1, 1]
void normalizeSignal(float* signal, size_t size) {
    float min_val = signal[0];
    float max_val = signal[0];
    
    for (size_t i = 1; i < size; i++) {
        if (signal[i] < min_val) min_val = signal[i];
        if (signal[i] > max_val) max_val = signal[i];
    }
    
    float range = max_val - min_val;
    if (range == 0) return;  // Avoid division by zero
    
    float scale = 2.0f / range;
    for (size_t i = 0; i < size; i++) {
        signal[i] = (signal[i] - min_val) * scale - 1.0f;
    }
}

// Simple moving average filter (window_size must be <= 16 in this snippet)
void movingAverageFilter(float* signal, size_t size, int window_size) {
    SimpleBuffer<float, 16> avgBuffer;
    float sum = 0;

    for (size_t i = 0; i < size; i++) {
        avgBuffer.push(signal[i]);
        sum = 0;
        for (size_t j = 0; j < avgBuffer.size(); j++) {
            sum += avgBuffer.get(j);
        }
        signal[i] = sum / avgBuffer.size();
    }
}

// Find top N=3 frequency peaks in vReal array (already in magnitude)
void findTop3Peaks(const float* magnitudeArray, size_t n, float& f1, float& f2, float& f3,
                   float& m1, float& m2, float& m3) 
{
    // For storing {index, magnitude}
    struct PeakInfo {
      size_t index;
      float mag;
    } topPeaks[TOP_PEAKS_COUNT] = {{0,0},{0,0},{0,0}};
    
    // Find top 3
    for (size_t i = 1; i < n; i++) {
        float currentMag = magnitudeArray[i];
        // Insert into topPeaks array if it is large enough
        for (int p = 0; p < TOP_PEAKS_COUNT; p++) {
            if (currentMag > topPeaks[p].mag) {
                // Shift lower peaks down
                for (int q = TOP_PEAKS_COUNT - 1; q > p; q--) {
                    topPeaks[q] = topPeaks[q-1];
                }
                // Insert current
                topPeaks[p].index = i;
                topPeaks[p].mag   = currentMag;
                break;
            }
        }
    }
    
    // Map to frequencies & dB magnitudes
    f1 = topPeaks[0].index * (float(SAMPLE_RATE) / (2.0f * float(n)));
    f2 = topPeaks[1].index * (float(SAMPLE_RATE) / (2.0f * float(n)));
    f3 = topPeaks[2].index * (float(SAMPLE_RATE) / (2.0f * float(n)));
    
    m1 = 20.0f * log10f(topPeaks[0].mag);
    m2 = 20.0f * log10f(topPeaks[1].mag);
    m3 = 20.0f * log10f(topPeaks[2].mag);
}

// Compute approximate PSD in dB over the entire signal
// PSD ~ 10 * log10( mean( x^2 ) )
float computePSD(const float* signal, size_t size) {
    double sumSq = 0;
    for (size_t i = 0; i < size; i++) {
        sumSq += (double)signal[i] * (double)signal[i];
    }
    float meanSq = (float)(sumSq / (double)size);
    return 10.0f * log10f(meanSq + 1.0e-12f); // add small epsilon to avoid log(0)
}

// Compute RMS of the signal
float computeRMS(const float* signal, size_t size) {
    double sumSq = 0;
    for (size_t i = 0; i < size; i++) {
        sumSq += (double)signal[i] * (double)signal[i];
    }
    float meanSq = (float)(sumSq / (double)size);
    return sqrtf(meanSq);
}

// Compute zero-crossing rate (ZCR)
float computeZCR(const float* signal, size_t size) {
    if (size < 2) return 0.0f;
    unsigned int count = 0;
    for (size_t i = 1; i < size; i++) {
        if (signOf(signal[i]) != signOf(signal[i-1])) {
            count++;
        }
    }
    return (float)count / (float)size;
}

void compressSoundFeatures(const struct SoundFeatures *input, struct SoundFeaturesCompressed *output) {
    output->dominantFrequency = (uint16_t)(input->dominantFrequency * 10);  // Scale factor 10
    output->magnitude = (uint16_t)(input->magnitude * 100);  // Scale factor 100
    output->frequencyDeviation = (uint16_t)(input->frequencyDeviation * 10);
    // output->timestamp = input->timestamp;
    
    output->topFrequency1 = (uint16_t)(input->topFrequency1 * 10);
    output->topFrequency2 = (uint16_t)(input->topFrequency2 * 10);
    output->topFrequency3 = (uint16_t)(input->topFrequency3 * 10);

    output->topMagnitude1 = (uint8_t)(input->topMagnitude1 * 10);
    output->topMagnitude2 = (uint8_t)(input->topMagnitude2 * 10);
    output->topMagnitude3 = (uint8_t)(input->topMagnitude3 * 10);
    
    output->anomalyCount = input->anomalyCount;  // No scaling needed
}

// ------------------ BLE Transmission -----------------------------
void transmitFeatures(const SoundFeatures& features) {
    if (deviceConnected) {
        struct SoundFeaturesCompressed compressed;
        compressSoundFeatures(&features, &compressed);
        pCharacteristic->setValue((uint8_t*)&compressed, sizeof(SoundFeaturesCompressed));
        pCharacteristic->notify();
    }
    else {
        Serial.println("Device not connected");
    }
}

// ------------------ Main Processing Function ----------------------
void processAndTransmitSignal() {
    // We'll copy the raw (original) samples into one array 
    // and the "processed" samples into processed_signal
    static float original_signal[SIGNAL_SIZE];

    // 1. Read from mock signal 
    for (size_t i = 0; i < SIGNAL_SIZE; i++) {
        float sampleVal = (float)mock_signal[mockSignalIndex];
        original_signal[i] = sampleVal;    // keep original
        processed_signal[i] = sampleVal;   // for further processing
        mockSignalIndex = (mockSignalIndex + 1) % SIGNAL_SIZE;  // wrap around
    }

    // 2. Normalize
    normalizeSignal(processed_signal, SIGNAL_SIZE);
    
    // 3. Apply moving average filter
    movingAverageFilter(processed_signal, SIGNAL_SIZE, 5);

    // 4. Prepare FFT
    for (size_t i = 0; i < SIGNAL_SIZE; i++) {
        vReal[i] = processed_signal[i];
        vImag[i] = 0;
    }

    // 5. Perform FFT
    //    (We use half of the bins for magnitude; the rest are mirrored)
    FFT.windowing(vReal, SIGNAL_SIZE, FFT_WIN_TYP_HAMMING, FFT_FORWARD);
    FFT.compute(vReal, vImag, SIGNAL_SIZE, FFT_FORWARD);
    FFT.complexToMagnitude(vReal, vImag, SIGNAL_SIZE);

    // 6. Find dominant frequency & magnitude in dB
    float maxMagnitude = 0;
    float dominantFrequency = 0;
    for (size_t i = 1; i < SIGNAL_SIZE / 2; i++) {
        if (vReal[i] > maxMagnitude) {
            maxMagnitude = vReal[i];
            dominantFrequency = i * ((float)SAMPLE_RATE / (float)SIGNAL_SIZE);
        }
    }
    float magnitudeDB = 20.0f * log10f(maxMagnitude + 1e-12f);

    // 7. Frequency deviation from 400Hz
    float freqDeviation = fabsf(dominantFrequency - TARGET_FREQUENCY);

    // 8. Detect anomaly if conditions are met 
    bool isAnomaly = false;
    if (freqDeviation > 50 || magnitudeDB > THRESHOLD_DB) {
        // Throttle anomaly counting to once per second
        if (millis() - stats.lastAnomalyTime > 1000) {
            isAnomaly = true;
            stats.anomalyCount++;
            stats.lastAnomalyTime = millis();
        }
    }

    // 9. Compute additional features (top 3 peaks, PSD, RMS, ZCR)
    float topF1=0, topF2=0, topF3=0;
    float topM1=0, topM2=0, topM3=0;
    findTop3Peaks(vReal, SIGNAL_SIZE/2, topF1, topF2, topF3, topM1, topM2, topM3);

    float psdValue = computePSD(original_signal, SIGNAL_SIZE);
    float rmsValue = computeRMS(original_signal, SIGNAL_SIZE);
    float zcrValue = computeZCR(original_signal, SIGNAL_SIZE);

    // 10. Prepare and transmit features
    SoundFeatures features;
    features.dominantFrequency  = dominantFrequency;
    features.magnitude          = magnitudeDB;
    features.frequencyDeviation = freqDeviation;
    features.anomalyCount       = stats.anomalyCount;
    // features.timestamp          = millis();

    features.topFrequency1  = topF1;
    features.topFrequency2  = topF2;
    features.topFrequency3  = topF3;
    features.topMagnitude1  = topM1;
    features.topMagnitude2  = topM2;
    features.topMagnitude3  = topM3;
    features.psd            = psdValue;
    features.rmsValue       = rmsValue;
    features.zcr            = zcrValue;

    transmitFeatures(features);

    // Debug output
    Serial.printf("Freq: %.1f Hz, Mag: %.1f dB, Dev: %.1f Hz%s\n",
                  dominantFrequency, magnitudeDB, freqDeviation,
                  isAnomaly ? " ⚠️ ANOMALY" : "");

    // Additional debug prints
    Serial.printf("Top peaks: (%.1f Hz, %.1f dB), (%.1f Hz, %.1f dB), (%.1f Hz, %.1f dB)\n",
                  topF1, topM1, topF2, topM2, topF3, topM3);
    Serial.printf("PSD: %.2f dB, RMS: %.3f, ZCR: %.3f\n\n",
                  psdValue, rmsValue, zcrValue);
}

// ------------------ BLE Setup and Arduino Setup ------------------
void setupBLE() {
    BLEDevice::init("InHive_Sensor");
    pServer = BLEDevice::createServer();
    pServer->setCallbacks(new ServerCallbacks());

    BLEService *pService = pServer->createService(SERVICE_UUID);
    pCharacteristic = pService->createCharacteristic(
        CHARACTERISTIC_UUID,
        BLECharacteristic::PROPERTY_NOTIFY | BLECharacteristic::PROPERTY_READ
    );
    
    pCharacteristic->addDescriptor(new BLE2902());
    pService->start();
    pServer->getAdvertising()->start();
    
    Serial.println("BLE Server ready for gateway connection");
}

void setup() {
    Serial.begin(115200);
    delay(1000);
    
    // Allocate memory for signal processing
    processed_signal = (float*)malloc(SIGNAL_SIZE * sizeof(float));
    vReal = (float*)malloc(SIGNAL_SIZE * sizeof(float));
    vImag = (float*)malloc(SIGNAL_SIZE * sizeof(float));
    
    if (!processed_signal || !vReal || !vImag) {
        Serial.println("Memory allocation failed!");
        return;
    }

    setupBLE();
    Serial.println("InHive Sensor Started - Using Mock Signal");
}

void loop() {
    if (millis() - lastUpdateTime >= UPDATE_INTERVAL) {
        processAndTransmitSignal();
        lastUpdateTime = millis();
    }
}