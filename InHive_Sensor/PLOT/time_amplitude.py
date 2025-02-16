import numpy as np
import matplotlib.pyplot as plt

# Load the mock signal data
signal = np.load("X:\Signal\InHive_Sensor\scripts\mock_signal.npy")

# Generate time axis assuming a uniform sampling rate
ts = 1 / 8000  # Example: 44.1 kHz sampling rate (modify if needed)
time = np.arange(len(signal)) * ts

# Plot the signal
plt.figure(figsize=(10, 4))
plt.plot(time, signal, label="Mock Signal")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.title("Continuous Mock Signal Plot")
plt.legend()
plt.grid()
plt.show()
