
# Repository Structure

- **`EE_vs_AMP.py`**: This script analyzes the energy efficiency of the decoder with respect to the number of AMP (Approximate Message Passing) iterations for different Signal-to-Noise Ratio (SNR) values.
- **`EE_vs_SNR.py`**:  This script analyzes the energy efficiency of the decoder against SNR for the optimal parameters (in terms of AMP and BP iterations) found in other simulations.
- **`EE_Decoder_vs_BlockLength.py`**: This script simulates the energy efficiency of the decoder with varying block lengths for different SNRs.
- **`EE_Encoder_vs_Blocklength.py`**: This script simulates the energy efficiency of the encoder with varying block lengths for different SNRs.
- **`generate_message.py`**: This script generates a random sequence of message bits, representing a message to be encoded and transmitted.
- **`script.sh`**: This is just to generate required number of messages to test the decoder/encoder
- **`srldpc.py`**: This script contains the core simulation functions for SPARC code encoding, decoding, and performance evaluation. 
- **`ldpc_graphs`**: This directory can hold LDPC code definitions in the form of `.alist` files. These files define the structure of the LDPC code used in the simulations. 
- **`messages`**: This directory holds `.txt` files containing sequences of message bits that are used as input to the encoder.
- **`simulations`**: This directory holds `.txt` files and plots of the simulations.

# Dependencies

- Python 3.7 or later
- NumPy
- pyfht
- gfldpc 
- matplotlib (for plotting)

You can install the necessary packages by running:

```bash
pip install -r requirements.txt
```

# Running the Simulations


## Generate Message Bits (Optional):

If you need a new random message sequence, run generate_message.py:

bash script.sh <num_messages>


This will create .txt files in the messages directory containing the message bits.

## Select and Execute a Simulation Script:

Choose the script corresponding to the experiment you want to run (e.g., EE_vs_AMP.py, EE_vs_SNR.py).

Modify simulation parameters within the script if needed (e.g., block_length, snr_value).

## Run the script:

python3 <selected_script>.py

Replace <selected_script>.py with the actual filename. You can also pass command-line arguments to override default parameters.

## Analyze Results:

Simulation results will be saved to files in the simulations directory.

Some scripts may generate plots, which will also be saved in the simulations directory.

# References 
https://github.com/EngProjects/mMTC.git