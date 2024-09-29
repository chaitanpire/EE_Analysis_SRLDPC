import os
import numpy as np
import srldpc
import argparse
from os import system

# Function to extract non-zero y_decode values and find the minimum value along with its index
def extract_min_y_decode(file_path):
    y_decode = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()  # Remove leading and trailing whitespace
            if line:  # Check if the line is not empty
                try:
                    y_decode.append(float(line))  # Try to convert to float
                except ValueError:
                    # Log the issue and skip the line if it is non-numeric
                    print(f"Skipping non-numeric line: {line}")

    # Convert to numpy array
    y_decode = np.array(y_decode)
    
    # Ignore zero values
    non_zero_y_decode = y_decode[y_decode > 0]
    
    # Return the minimum value and its index if non-zero values exist, else None
    if len(non_zero_y_decode) > 0:
        min_value = np.min(non_zero_y_decode)
        min_index = np.where(y_decode == min_value)[0][0]  # Get index of min value
        return min_value, min_index
    else:
        return None, None

# Directory where all your files are located
directory = "./simulations"  # Update this to the correct directory path

# Dictionary to store minimum y_decode, BP, and AMP for each SNR
min_y_decode_info_per_snr = {}

# Iterate through all files in the directory
for file_name in os.listdir(directory):
    if file_name.startswith("snr_") and file_name.endswith("_y_decode.txt"):
        # Extract SNR and BP values from the file name
        parts = file_name.split('_')
        snr_value = float(parts[1])
        bp_value = float(parts[3])
        
        # Construct the full file path
        file_path = os.path.join(directory, file_name)
        
        # Extract the minimum y_decode value and its index (AMP)
        min_y_decode, min_index = extract_min_y_decode(file_path)
        
        # AMP is inferred from the line number (index)
        if min_y_decode is not None:
            amp_value = min_index  # AMP is inferred from the line number (index)
            
            # Store the minimum y_decode, BP, AMP, and index for the corresponding SNR
            if snr_value not in min_y_decode_info_per_snr:
                min_y_decode_info_per_snr[snr_value] = {
                    'min_y_decode': min_y_decode,
                    'bp': bp_value,
                    'amp': amp_value,
                    'file': file_name
                }
            else:
                # Update only if a new lower y_decode is found for the SNR
                if min_y_decode < min_y_decode_info_per_snr[snr_value]['min_y_decode']:
                    min_y_decode_info_per_snr[snr_value] = {
                        'min_y_decode': min_y_decode,
                        'bp': bp_value,
                        'amp': amp_value,
                        'file': file_name
                    }

# Print the extracted values for each SNR
print("Minimum y_decode values along with BP and AMP for each SNR:")
for snr, info in min_y_decode_info_per_snr.items():
    print(f"SNR: {snr}, Min y_decode: {info['min_y_decode']}, BP: {info['bp']}, AMP: {info['amp']}, File: {info['file']}")
#  Now run the simulations for the given SNR and BP and AMP values for 


# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--keep_graph', type=int, default=1,
                    help='Flag for keeping factor graph messages between AMP iterations')
parser.add_argument('--snr_value', type=int, default=3)
parser.add_argument('--block_length', type=int, default=7350)

args = parser.parse_args()

# Extract CL arguments
max_num_amp_iter = 50
num_final_bp_iter = 300
keep_graph = args.keep_graph
snr_value = args.snr_value
block_length = args.block_length
num_chnl_uses = 7350

# Define simulation directory
sim_dir = 'simulations'
system(f'mkdir -p ./{sim_dir}')

# Select Eb/N0 range
ebno_db_vals = np.arange(2.0, 3.01, 0.25)


# Execute simulation
# Loop through the messages and perform simulation
# Loop through the messages and perform simulation
# Loop through the messages and perform simulation
list_messages = os.listdir("./messages/")
num_iter = 1
log = 0

# Get the number of SNR values from the dictionary
num_snr_values = len(min_y_decode_info_per_snr)
y_decode_final = np.zeros(num_snr_values)

# Extract SNR values from the dictionary and sort them
sorted_snr_info = sorted(min_y_decode_info_per_snr.items())  # Sort by SNR

# Prepare arrays for SNR and y_decode values
snr_values_sorted = np.array([float(snr) for snr, info in sorted_snr_info])
y_decode_sorted = np.zeros_like(snr_values_sorted)

# Directory for storing simulation results
sim_dir2 = 'results/revisedbler'

# Loop through the messages and iterate over each sorted SNR value
for messages in list_messages:
    m_file = f"./messages/{messages}"
    for iter in range(num_iter):
        z = np.random.randn(num_chnl_uses, 1)
        for i, (snr, info) in enumerate(sorted_snr_info):  # Use sorted SNR values
            amp_value = int(info['amp'])
            bp_value = int(info['bp'])
            
            # Ensure snr is a float (already handled in sorting)
            snr_float = float(snr)
            
            print(f"Running simulation for SNR: {snr}, BP: {bp_value}, AMP: {amp_value}, message file: {messages}, iteration: {iter}")
            # Pass the float SNR to the simulation
            x_, y_ = srldpc.simulate(
                z, m_file, snr_float, sim_dir2, log, num_chnl_uses, amp_value,
                num_final_bp_iter, bp_value, keep_graph
            )
            y_decode_sorted[i] += y_

# Average the results over the number of iterations and messages
y_decode_sorted /= (num_iter * len(list_messages))

# Now plot these sorted y_decode values against the sorted SNR values
import matplotlib.pyplot as plt

plt.plot(snr_values_sorted, 5888 / (28 * y_decode_sorted), marker = 'o')  # The formula for energy efficiency
plt.title('Energy Efficiency Analysis of the Decoder')
plt.xlabel('SNR')
plt.ylabel('Energy Efficiency (bits/J)')

plt.savefig('./simulations/plots/EE_vs_SNR_plot.png')
