"""
Sparse Regression LDPC simulation.
"""

import numpy as np
from time import time
from os import system
from os.path import join
from pyfht_edit import block_sub_fht
from gfldpc import GFLDPC
import matplotlib as plt

__author__ = 'Ebert'
__date__ = '18 January 2023'

def sparc_codebook(num_chnl_uses, q, n_gf):
    """
    Create functions for efficiently computing Ax and A^Tz. 

        Arguments:
            num_chnl_uses (int): number of channel uses (real d.o.f.)
            q (int): length of each section
            n_gf (int): number of sections

        Returns:
            Ab(b) (function): computes matrix-vector product Ab
            Az(z) (function): computes matrix-vector product A^Tz
    """
    Ax, Ay, _ = block_sub_fht(num_chnl_uses, q, n_gf, seed=None, ordering=None)
    def Ab(b):
        return Ax(b).reshape(-1, 1)/np.sqrt(num_chnl_uses)
    def Az(z):
        return Ay(z).reshape(-1, 1)/np.sqrt(num_chnl_uses) 
    return Ab, Az

def compute_likelihood_vector(r0, tau, d, n, q):
    """
    Compute likelihood vector according to eq (10). 

        Arguments:
            r0 (LMx1 ndarray): AMP effective observation
            tau (float): AMP effective noise variance parameter
            d (float): amplitude scaling
            n (int): number of sections
            q (int): length of each section

        Returns:
            alpha (MxL ndarray): likelihood vectors for all L sections

    """
    r = r0.copy().reshape(n, q)
    alpha = np.exp(d*r/tau**2)
    row_norms = np.linalg.norm(alpha, ord=1, axis=1).reshape(-1, 1)
    return alpha / row_norms

def bp_denoiser(r, code, tau, d, keep_graph, num_bp_iter=1):
    """
    BP Denoiser from Sparse Regression LDPC Codes.

        Arguments:
            r (LMx1 ndarray): AMP effective observation
            code (object): LDPC outer code
            tau (float): AMP effective noise variance parameter
            d (float): amplitude scaling
            keep_graph (bool): flag of whether to retain graph messages
            num_bp_iter (int): number of BP iterations to perform

        Returns:
            q (LMx1 ndarray): denoised state estimate
    """
    n = code.N
    q = code.q

    if keep_graph:
        code.reset_observations()
    else:
        code.reset_graph()

    alpha = compute_likelihood_vector(r, tau, d, n, q)
    for i in range(n):
        code.set_observation(i, alpha[i, :])
    code.bp_decoder(num_bp_iter)
    
    return code.get_estimates().reshape(-1, 1)

def amp_state_update(z, s, d, Az, num_bp_iter, keep_graph, code):
    """
    Compute state update within AMP iterate. 

        Arguments:
            z (nx1 ndarray): AMP residual vector
            s (LMx1 ndarray): AMP state estimate
            d (float): amplitude scaling per entry
            Az (function): SPARC encoding function
            num_bp_iter (int): number of BP iterations to perform
            keep_graph (bool): flag of whether to retain graph messages
            code (object): outer LDPC code object

        Returns:
            s_plus (LMx1 ndarray): updated state estimate
    """
    n = z.size
    tau = np.sqrt(np.sum(z**2)/n)
    r = (d*s + Az(z))
    s_plus = bp_denoiser(r, code, tau, d, keep_graph, num_bp_iter)

    return s_plus

def amp_residual(y, z, s, d, Ab):
    """
    Compute residual within AMP iterate.

        Arguments: 
            y (nx1 ndarray): received vector
            z (nx1 ndarray): AMP residual vector
            s (LMx1 ndarray): current AMP state estimate
            d (float): amplitude scaling per entry
            Ab (function): SPARC encoding function

        Returns:
            z (nx1 ndarray): updated AMP residual vector
    """
    n = y.size
    tau = np.sqrt(np.sum(z**2)/n)
    onsager_term = (d**2)*(np.sum(s) - np.sum(s**2))
    z_plus = y - d*Ab(s) + (z/(n*tau**2))*onsager_term
    
    return z_plus
def quantise( x, q):
        """
        Quantise the input signal x to q levels. 

            Arguments:
                x (nx1 ndarray): input signal
                q (int): number of quantisation levels

            Returns:
                xq (nx1 ndarray): quantised signal
        """
        xq = np.zeros(x.shape)
        for i in range(x.size):
            xq[i] = np.round((q-1)*x[i])/q-1
        return xq
def simulate(z,message_file,ebnodb, 
             sim_dir, loglist: list,num_chnl_uses,
             num_amp_iter, 
             num_final_bp_iter, 
             num_bp_denoiser_iter, 
             keep_graph=True):
    
    # Define LDPC code over finite field by specifying the alist filename
    code_folder = './ldpc_graphs/varying_ldpc_rates/'
    code_filename = 'ldpc_n766_k736_gf256.txt'
    code = GFLDPC(join(code_folder, code_filename))

    # Define filename
    filename = f'./{sim_dir}/ebnodb_{ebnodb}_2_numampiter_{num_amp_iter}_' + \
               f'numfinalbpiter_{num_final_bp_iter}_numbpdenoiseriter_' + \
               f'{num_bp_denoiser_iter}_keepgraph_{keep_graph}_' + \
               f'{code_filename[:-4]}.txt'

    # Code parameters
    q = code.q
    bits_per_symbol = int(np.log2(q))
    k_gf = code.K
    n_gf = code.N
    k_bin = k_gf*bits_per_symbol
    # num_chnl_uses = 7350

    # Decoder parameters
    target_num_block_errors = 50

    # Channel noise parameters
    ebno = 10**(ebnodb/10)
    nvar, sigma = 1, 1
    noise_psd = 2 * nvar
    energy_per_bit = ebno * noise_psd
    total_energy = energy_per_bit * k_bin
    column_energy_scaling = total_energy / n_gf
    d = np.sqrt(column_energy_scaling)

    # Printing options
    should_print = True
    write_frequency = 50

    # Error-tracking data structures
    ber = 0.0
    bler = 0.0
    num_bit_errors = 0
    num_block_errors = 0
    num_sims = 0
    snr = 0.0


    # Monte-carlo simulation to estimate BER/BLER
    while (num_block_errors < target_num_block_errors):
        
        # Periodically update data
        if (should_print and (num_sims % write_frequency == 0)):
            log = open(filename, 'w')
            log.write(f'BER:\t{ber:.4e}\n')
            log.write(f'BLER:\t{bler:.4e}\n')
            log.write(f'Num Bit Errors:\t{num_bit_errors}\n')
            log.write(f'Num Block Errors:\t{num_block_errors}\n')
            log.write(f'Completed sims:\t{num_sims}\n')
            log.write(f'Empirical SNR:\t{snr}\n')
            log.close()
        
        # Reset SRLDPC outer factor graph
        code.reset_graph()

        # Randomly select LDPC codeword
        # bin_msg = np.random.randint(2, size=k_bin)
        #read message from file message.txt
        with open(message_file, 'r') as file:
            bin_msg = file.read()
        bin_msg = bin_msg.split()
        
        bin_msg = [int(i) for i in bin_msg]
        bin_msg = np.array(bin_msg)
        assert bin_msg.size == k_bin
        # Encode message
        tic_encode = time()

        codeword = code.encode(bin_msg)
        assert code.check_consistency(codeword)

        # Generate sparse representation through indexing
        sparse_codeword = np.zeros(q*n_gf)
        sparse_codeword[np.arange(n_gf)*q+codeword] = 1        
        

        # Generate the binned SPARC codebook
        Ab, Az = sparc_codebook(num_chnl_uses, q, n_gf)
        # Generate transmitted signal
        x = d*Ab(sparse_codeword)
        quantised = quantise(x, 65536)

        toc_encode = time()
        
        # Send signal through AWGN channel
        z = sigma*z
        z = z.reshape(-1, 1)
        
        y = x + z
        nvarht = np.linalg.norm(z)**2 / num_chnl_uses
        snr = (num_sims*snr + 10*np.log10(np.linalg.norm(x)**2 / 
                                          (k_bin*2*nvarht))) / (num_sims+1)

        tic_decode = time()
        # Prepare for AMP decoding
        z = y.copy()
        s = np.zeros((q*n_gf, 1))
        # AMP decoding
        for idx_amp_iter in range(num_amp_iter):
            # Adjust number of BP iterations
            if num_bp_denoiser_iter == -1: 
                num_bp_iter = idx_amp_iter + 1
            else:
                num_bp_iter = num_bp_denoiser_iter

            # AMP iterate
            s = amp_state_update(z, s, d, Az, num_bp_iter, keep_graph, code)
            z = amp_residual(y, z, s, d, Ab)
            
            # Check stopping conditions
            cdwd_ht = np.array([np.argmax(s[i*q:(i+1)*q]) for i in range(n_gf)])
            if code.check_consistency(cdwd_ht):
                break  
            # Run BP decoder
            if idx_amp_iter == (num_amp_iter - 1):
                s = amp_state_update(z, s, d, Az, num_final_bp_iter, \
                                     False, code)

        # Make hard decisions on final estimated vector 
        s = s.reshape(n_gf, q)
        codeword_ht = np.argmax(s, axis=1).flatten().astype(int)
        toc_decode = time()
        # Compute BER/BLER  
        update_bler = True
        for i in range(n_gf):
            err_i = codeword[i]^codeword_ht[i]
            binstr = '0'+ bin(err_i)[2:]
            err_i_bits = [int(bit) for bit in binstr]
            if i >= (n_gf - k_gf):
                num_bit_errors += np.sum(err_i_bits)
            if update_bler and (np.sum(err_i_bits) > 0):
                num_block_errors += 1
                update_bler = False
        
        # Increment error-tracking parameters
        num_sims += 1
        ber = num_bit_errors/(num_sims*k_bin)
        bler = num_block_errors/(num_sims)
        print("Num_BP_Iter: ", num_bp_denoiser_iter)    
        print("Decode Time: ", toc_decode - tic_decode)

    # Record final results.
    print(f'COMPLETED Eb/No: {ebnodb}dB\tBER: {ber}\tBLER: {bler}')
    log = open(filename, 'w')
    log.write(f'BER:\t{ber:.4e}\n')
    log.write(f'BLER:\t{bler:.4e}\n')
    log.write(f'Num Bit Errors:\t{num_bit_errors}\n')
    log.write(f'Num Block Errors:\t{num_block_errors}\n')
    log.write(f'Completed sims:\t{num_sims}\n')
    log.write(f'Empirical SNR:\t{snr}\n')
    log.close()
    return toc_encode - tic_encode, toc_decode - tic_decode 

