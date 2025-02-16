import numpy as np
import os

# Define the correct path for mock_signal.npy
script_dir = os.path.dirname(__file__)  # Get the directory where the script is located
npy_path = os.path.join(script_dir, "X:\Signal\InHive_Sensor\scripts\mock_signal.npy")  # Full path to the .npy file

# Check if the file exists
if not os.path.exists(npy_path):
    print(f"Error: {npy_path} not found!")
    exit(1)

# Load the mock signal
mock_signal = np.load(npy_path)

# Offset the signal so that all values are positive (prevent negative clipping)
offset = np.abs(np.min(mock_signal))  # Find the lowest negative value and shift the signal
mock_signal_shifted = mock_signal + offset  # Shift so the minimum is 0

# Scale up to fit uint16 range (to preserve precision without normalizing)
mock_signal_int = (mock_signal_shifted * 4095).astype(np.uint16)

# Save as C header inside the include/ directory
header_path = os.path.join(script_dir, "../include/mock_signal.h")
with open(header_path, "w") as f:
    f.write("#ifndef MOCK_SIGNAL_H\n#define MOCK_SIGNAL_H\n\n")
    f.write("const uint16_t mock_signal[] = {")
    f.write(",".join(map(str, mock_signal_int)))
    f.write("};\n")
    f.write(f"const int mock_signal_size = {len(mock_signal_int)};\n")
    f.write("#endif")

print(f"mock_signal.h successfully saved to {header_path}")
