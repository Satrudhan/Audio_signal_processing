// InHive_Gateway.ino

#include <Arduino.h>
#include <BLEDevice.h>
#include <BLEUtils.h>
#include <BLEScan.h>
#include <BLEAdvertisedDevice.h>
#include <WiFi.h>
#include <HTTPClient.h>
#include <ArduinoJson.h>

// ---------------- WiFi credentials ----------------
const char* ssid = "Sumanas";
const char* password = "Nepenthe108";

// ---------------- Server details ------------------
const char* serverUrl = "http://192.168.49.93:5000/receive_data";

// ---------------- BLE Constants -------------------
#define SERVICE_UUID        "91bad492-b950-4226-aa2b-4ede9fa42f59"
#define CHARACTERISTIC_UUID "cba1d466-344c-4be3-ab3f-189f80dd7518"

// -------------------------------------------------------------------------------
// Updated data structure matching the new sensor firmware
// Note: Must match exactly in size and field order so memcpy works correctly.
// -------------------------------------------------------------------------------
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

// ---------------- Global variables for BLE -----------------------
BLEClient* pClient = nullptr;
BLERemoteCharacteristic* pRemoteCharacteristic = nullptr;

// ---------------- Global connection status tracking --------------
struct ConnectionStatus {
  bool bleConnected = false;
  bool wifiConnected = false;
  unsigned long lastDataReceived = 0;
  unsigned long lastServerSync = 0;
  int failedUploads = 0;
} status;

// ---------------- Queue for storing data when offline ------------
const int QUEUE_SIZE = 50;
SoundFeatures dataQueue[QUEUE_SIZE];
int queueHead = 0;
int queueSize = 0;

// WiFi check every 30 seconds
unsigned long lastWiFiCheck = 0;
const int WiFiCheckInterval = 30000;

// Forward declaration for sendToServer
bool sendToServer(const SoundFeatures& features);

// ---------------- BLE Client Callbacks ---------------------------
class ClientCallbacks : public BLEClientCallbacks {
  void onConnect(BLEClient* pclient) override {
    status.bleConnected = true;
    Serial.println("Connected to InHive sensor!");
  }
  
  void onDisconnect(BLEClient* pclient) override {
    status.bleConnected = false;
    Serial.println("Disconnected from InHive sensor!");
  }
};

void decompressSoundFeatures(const struct SoundFeaturesCompressed *input, struct SoundFeatures *output) {
  output->dominantFrequency = (float)input->dominantFrequency / 10.0;
  output->magnitude = (float)input->magnitude / 100.0;
  output->frequencyDeviation = (float)input->frequencyDeviation / 10.0;

  output->topFrequency1 = (float)input->topFrequency1 / 10.0;
  output->topFrequency2 = (float)input->topFrequency2 / 10.0;
  output->topFrequency3 = (float)input->topFrequency3 / 10.0;

  output->topMagnitude1 = (float)input->topMagnitude1 / 10.0;
  output->topMagnitude2 = (float)input->topMagnitude2 / 10.0;
  output->topMagnitude3 = (float)input->topMagnitude3 / 10.0;

  output->anomalyCount = input->anomalyCount;
}

// ---------------- BLE Notification Callback ----------------------
// Called every time a BLE notification is received.
void notifyCallback(BLERemoteCharacteristic* pBLERemoteCharacteristic, 
                    uint8_t* pData, size_t length, bool isNotify) 
{
  // Ensure the received data length matches our struct
  Serial.println(sizeof(SoundFeaturesCompressed));
  Serial.println(length);
  Serial.println(sizeof(pData));

  if (length == sizeof(SoundFeaturesCompressed)) {
    SoundFeaturesCompressed features;
    SoundFeatures decompressedFeatures;
    memcpy(&features, pData, sizeof(SoundFeaturesCompressed));
    status.lastDataReceived = millis();
    decompressSoundFeatures(&features, &decompressedFeatures);
    // Print a quick debug message
    Serial.println("Notification received from sensor!");
    
    // Attempt sending it directly to the server
    if (!sendToServer(decompressedFeatures)) {
      // If the send fails, add to the circular queue (if space exists)
      if (queueSize < QUEUE_SIZE) {
        dataQueue[(queueHead + queueSize) % QUEUE_SIZE] = decompressedFeatures;
        queueSize++;
        Serial.println("Data queued for later transmission.");
      } else {
        // If queue is full, overwrite the oldest data
        Serial.println("Data queue full, dropping oldest entry.");
        queueHead = (queueHead + 1) % QUEUE_SIZE;
        dataQueue[(queueHead + queueSize - 1) % QUEUE_SIZE] = decompressedFeatures;
      }
    }
  }
  else {
    Serial.println("length mismatch");
  }
}

// ---------------- Process Queue (Resend Data) --------------------
void processQueue() {
  while (queueSize > 0 && WiFi.status() == WL_CONNECTED) {
    if (sendToServer(dataQueue[queueHead])) {
      // Successfully sent. Remove from the queue.
      queueHead = (queueHead + 1) % QUEUE_SIZE;
      queueSize--;
    } else {
      // If a send fails, exit processing for now, we'll retry later.
      break;
    }
  }
}

// ---------------- Sends the sensor data to the server ------------
bool sendToServer(const SoundFeatures& features) {
  // Check if WiFi is connected; if not, bail early.
  if (WiFi.status() != WL_CONNECTED) {
    Serial.println("WiFi not connected. Cannot send data to server.");
    return false;
  }
  
  HTTPClient http;
  // Attempt connection to the server
  if (!http.begin(serverUrl)) {
    Serial.println("Failed to connect to server URL");
    return false;
  }
  http.addHeader("Content-Type", "application/json");
  
  // Create JSON document and serialize sensor data.
  // Bump this up from 256 to 512 to ensure enough space for additional fields.
  StaticJsonDocument<512> doc;
  doc["dominantFrequency"]   = features.dominantFrequency;
  doc["magnitude"]           = features.magnitude;
  doc["frequencyDeviation"]  = features.frequencyDeviation;
  doc["anomalyCount"]        = features.anomalyCount;
  // doc["timestamp"]           = features.timestamp;

  // Additional fields
  doc["topFrequency1"]   = features.topFrequency1;
  doc["topFrequency2"]   = features.topFrequency2;
  doc["topFrequency3"]   = features.topFrequency3;
  doc["topMagnitude1"]   = features.topMagnitude1;
  doc["topMagnitude2"]   = features.topMagnitude2;
  doc["topMagnitude3"]   = features.topMagnitude3;
  doc["psd"]             = features.psd;
  doc["rmsValue"]        = features.rmsValue;
  doc["zcr"]             = features.zcr;
  
  String jsonPayload;
  serializeJson(doc, jsonPayload);
  
  // POST the JSON to the server
  int httpResponseCode = http.POST(jsonPayload);
  if (httpResponseCode > 0) {
    // Successful POST, print result
    Serial.printf("POST successful! Code: %d\n", httpResponseCode);
    status.lastServerSync = millis();
    http.end();
    return true;
  } else {
    // Failed POST
    Serial.printf("POST failed. Error: %s\n", http.errorToString(httpResponseCode).c_str());
    status.failedUploads++;
    http.end();
    return false;
  }
}

// ---------------- Attempt to connect to BLE sensor ---------------
bool connectToBLE() {
  // Initialize BLE with a custom gateway name
  BLEDevice::init("InHive_Gateway");
  
  // Get and configure the BLE scanner
  BLEScan* pBLEScan = BLEDevice::getScan();
  pBLEScan->setActiveScan(true);
  Serial.println("Scanning for BLE devices...");
  
  // Perform a scan for 5 seconds
  BLEScanResults scanResults = pBLEScan->start(5);
  
  // Check if any devices were found
  if (scanResults.getCount() == 0) {
    Serial.println("No BLE devices found!");
    return false;
  }
  
  // List all detected devices and look for one named "InHive_Sensor"
  Serial.println("Devices found:");
  BLEAdvertisedDevice* myDevice = nullptr;
  for (int i = 0; i < scanResults.getCount(); i++) {
    BLEAdvertisedDevice device = scanResults.getDevice(i);
    
    // Get device name if available
    std::string deviceName = device.getName();
    String nameStr = deviceName.length() > 0 ? String(deviceName.c_str()) : "No Name";
    
    // Get service UUID if available (for logging)
    String serviceUUID = (device.haveServiceUUID()) 
                          ? String(device.getServiceUUID().toString().c_str()) 
                          : "None";
    
    Serial.printf("%d: Name: %s, Address: %s, RSSI: %d, Service UUID: %s\n",
                  i,
                  nameStr.c_str(),
                  device.getAddress().toString().c_str(),
                  device.getRSSI(),
                  serviceUUID.c_str());
    
    // Check if this device has the name "InHive_Sensor"
    if (nameStr == "InHive_Sensor") {
      myDevice = new BLEAdvertisedDevice(device);
      // If you want to break on the first match:
      // break;
    }
  }
  
  if (myDevice == nullptr) {
    Serial.println("InHive_Sensor not found by name");
    return false;
  }
  
  // Create a BLE client and set callbacks.
  pClient = BLEDevice::createClient();
  pClient->setClientCallbacks(new ClientCallbacks());
  
  Serial.println("Connecting to InHive_Sensor...");
  if (!pClient->connect(myDevice)) {
    Serial.println("BLE sensor connection failed");
    return false;
  }
  
  // Once connected, retrieve the service using the predefined SERVICE_UUID.
  BLERemoteService* pRemoteService = pClient->getService(SERVICE_UUID);
  if (pRemoteService == nullptr) {
    Serial.println("Failed to find the sensor service");
    pClient->disconnect();
    return false;
  }
  
  // Retrieve the characteristic using the predefined CHARACTERISTIC_UUID.
  pRemoteCharacteristic = pRemoteService->getCharacteristic(CHARACTERISTIC_UUID);
  if (pRemoteCharacteristic == nullptr) {
    Serial.println("Failed to find the sensor characteristic");
    pClient->disconnect();
    return false;
  }
  
  // Register for notifications if the characteristic supports it.
  if (pRemoteCharacteristic->canNotify()) {
    pRemoteCharacteristic->registerForNotify(notifyCallback);
    Serial.println("Notifications registered");
  } else {
    Serial.println("Characteristic does not support notifications");
    pClient->disconnect();
    return false;
  }
  
  return true;
}

// ---------------- Arduino Setup ----------------------------------
void setup() {
  Serial.begin(115200);
  delay(1000);
  
  // Connect to WiFi
  Serial.printf("Connecting to WiFi: %s\n", ssid);
  WiFi.begin(ssid, password);
  int wifiRetries = 0;
  while (WiFi.status() != WL_CONNECTED && wifiRetries < 20) {
    delay(500);
    Serial.print(".");
    wifiRetries++;
  }
  if (WiFi.status() == WL_CONNECTED) {
    Serial.println("\nConnected to WiFi!");
    status.wifiConnected = true;
  } else {
    Serial.println("\nWiFi connection failed. Will retry in loop.");
    status.wifiConnected = false;
  }
  
  // Connect to BLE sensor
  if (!connectToBLE()) {
    Serial.println("Initial BLE connection failed. Will attempt reconnect in loop.");
  }
}

// ---------------- Arduino Loop -----------------------------------
void loop() {
  unsigned long now = millis();
  
  // Check and maintain WiFi connection every WiFiCheckInterval ms
  if (now - lastWiFiCheck >= WiFiCheckInterval) {
    lastWiFiCheck = now;
    if (WiFi.status() != WL_CONNECTED) {
      Serial.println("WiFi lost. Attempting reconnection...");
      WiFi.disconnect();
      WiFi.begin(ssid, password);
      int attempts = 0;
      while (WiFi.status() != WL_CONNECTED && attempts < 20) {
        delay(500);
        Serial.print(".");
        attempts++;
      }
      if (WiFi.status() == WL_CONNECTED) {
        Serial.println("\nWiFi reconnected!");
        status.wifiConnected = true;
      } else {
        Serial.println("\nWiFi still disconnected.");
        status.wifiConnected = false;
      }
    }
  }
  
  // If BLE is disconnected, try to reconnect once in a while
  if (!status.bleConnected) {
    Serial.println("BLE sensor disconnected. Retrying connection...");
    if (connectToBLE()) {
      Serial.println("Reconnected to BLE sensor");
    } else {
      Serial.println("BLE reconnect failed. Will try again later.");
    }
    // Wait a bit before the next attempt
    delay(5000);
  }
  
  // Process queued sensor data if any
  if (status.wifiConnected && queueSize > 0) {
    processQueue();
  }
  
  // Small delay to avoid saturating CPU
  delay(10);
}