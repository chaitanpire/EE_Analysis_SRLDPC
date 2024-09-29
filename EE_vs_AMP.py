# EE_vs_AMP.py

import srldpc
import argparse
import numpy as np
from os import system
import os
import matplotlib.pyplot as plt

# Parse command line arguments
parser = argparse.ArgumentParser()


parser.add_argument('--num_bp_denoiser_iter', type=int, default=5,
                    help='Number of BP iter to perform per AMP iteration.')
parser.add_argument('--keep_graph', type=int, default=1,
                    help='Flag for keeping factor graph messages between AMP iterations')
parser.add_argument('--block_length', type=int, default=7350)

args = parser.parse_args()

# Extract CL arguments
max_num_amp_iter = 50
num_final_bp_iter = 300
keep_graph = args.keep_graph
block_length = args.block_length
num_chnl_uses = 7350

# Define simulation directory
sim_dir = 'simulations'
system(f'mkdir -p ./{sim_dir}')

# Select Eb/N0 range
ebno_db_vals = np.arange(2.0, 3.01, 0.25)

# Execute simulation
list_messages = os.listdir("./messages/")
max_num_bp_iter = 3
num_iter = 10
log = 0
sim_dir2 = 'results/revisedbler'
for snr in ebno_db_vals:
    y_decode = np.zeros([max_num_bp_iter+1, max_num_amp_iter + 1])  # Size based on max AMP iter
    x = np.arange(0, max_num_amp_iter + 1)
    for message_file in list_messages:
        m_file = f"./messages/{message_file}"
        # Initialize y_decode for the current BP iteration and AMP iterations
        for i in range(num_iter):
            z = np.random.randn(num_chnl_uses, 1)
            for num_bp_iter in range(0, max_num_bp_iter + 1):
                # Running simulation for each message file
                for num_amp_iter in range(5, max_num_amp_iter + 1, 5):
                    print(f"Running simulation for SNR: {snr}, BP iter: {num_bp_iter}, message file: {message_file}, iteration: {i}, AMP iteration: {num_amp_iter}")
                    # Run simulation and accumulate results based on AMP iterations
                    x_, y_ = srldpc.simulate(z, m_file, snr, sim_dir2, log, block_length,
                                                             num_amp_iter, num_final_bp_iter,
                                                             num_bp_iter, keep_graph)
                    y_decode[num_bp_iter][num_amp_iter] += y_
        # Average the results over the number of iterations
    y_decode /= (num_iter * len(list_messages))
    for num_bp_iter in range(0, max_num_bp_iter + 1):
        # Save results for each SNR and BP iter
        np.savetxt(f'./{sim_dir}/snr_{snr}_bp_{num_bp_iter}_y_decode.txt', y_decode[num_bp_iter])
        np.savetxt(f'./{sim_dir}/snr_{snr}_bp_{num_bp_iter}_x.txt', x)

    # Generate plots for each SNR
    for snr in ebno_db_vals:
        plt.figure()

        # Store y_decode for each num_amp_iter
        y_decode_list = []
        num_amp_iters = []

        for num_bp_iter in range(0, 4):  # Assuming max_num_bp_iter is 3
            y_decode_file = f'{sim_dir}/snr_{snr}_bp_{num_bp_iter}_y_decode.txt'
            if os.path.exists(y_decode_file):
                y_decode = np.loadtxt(y_decode_file)
                print(f"Loaded y_decode for SNR {snr}, BP {num_bp_iter}: {y_decode}")

                # Drop zero values
                non_zero_indices = y_decode > 0
                y_decode_filtered = y_decode[non_zero_indices]
                num_amp_iter_filtered = np.arange(0, len(y_decode))[non_zero_indices]

                print(f"Filtered y_decode for SNR {snr}, BP {num_bp_iter}: {y_decode_filtered}")
                print(f"The values have been plotted for AMP iterations: {num_amp_iter_filtered}")

                # Check if we have valid data to plot
                if len(y_decode_filtered) > 0:
                    y_decode_list.append(y_decode_filtered)
                    num_amp_iters.append(num_amp_iter_filtered)

        # Check if we have any valid data to plot
        if not y_decode_list:
            print(f"No valid data to plot for SNR {snr}. Skipping.")
            continue

        # Plot each y_decode
        for i, y_decode in enumerate(y_decode_list):
            plt.plot(num_amp_iters[i], 5888/(28*y_decode), label=f'BP Iteration: {i}')

        # Customize the plot
        plt.title(f'Energy Efficiency Analysis of the Decoder (SNR = {snr})')
        plt.xlabel('Number of AMP Iterations')
        plt.ylabel('Energy Efficiency (bits/J)')
        plt.legend()
        plt.grid()

        # Save the plot
        plt.savefig(f'{sim_dir}/plots/snr_{snr}_y_decode_vs_amp_iterations.png')
        plt.close()  # Close the plot to save memory