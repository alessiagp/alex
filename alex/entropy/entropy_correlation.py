import pandas as pd
import numpy as np
import scipy.integrate as integral
    
#------------------------------------------------------------------------------------------------#

def autocorrelation_function(f):
    """
    Function calculating autocorrelation for a given list, returns the autocorrelation function as a list
    """
    N = len(f) # total number of snapshots
    
    mean_f = np.mean(f)
    var_f = np.var(f)
    # k = t'/Dt
    
    Cf = []
    for k in range(N): # can cycle among different possible delays goes from 0 to N-1
        num = 0
        for j in range(N-k):
            num += (f[j] - mean_f) * (f[j+k] - mean_f)
        Cf_k = num/(N * var_f)
        Cf.append(Cf_k)
        
    return Cf
    
#------------------------------------------------------------------------------------------------#

def integrate_acf(autocorr, entropy_mean, mode="full"):
    """
    Function performing the integration of the autocorrelation function, with the option to choose whether to integrate it fully or just the positive part.
    Default is `full`.
    Returns a tuple with autocorrelation time, the number of independent samples, and the standard deviation of autocorrelation.
    """
    N = len(autocorr) # total number of elements in the autocorrelation function - that corresponds to the total number of snapshots
    h = 0 # threshold value
    ii = 0
    temp = 0
    eps = 1e-2
    time=range(0,N)
    t_sim = time[N-1]
    
    # Search for the autocorrelation time
    while (ii != N-1 and temp != 1):
        if((autocorr[ii] < (h+eps)) and (autocorr[ii] > (h-eps))):
            temp = 1
            index = ii
        else:
            ii = ii+1

    if mode == "positive":
        N_max = ii
        print(f"\nPOSITIVE ACF\nAutocorrelation function around 0: {autocorr[N_max]}\nN_max: {N_max}")
        t_auto_integ = integral.simps(autocorr[:N_max], time[:N_max])
        print(f"Autocorrelation time: {t_auto_integ} MC steps")

        N_ind_integ = t_sim/t_auto_integ
        SE_integ =  np.sqrt(np.var(entropy_mean)/N_ind_integ) # std deviation of the mean
        print(f"Number of independent samples: {N_ind_integ}\nStandard deviation of the mean: {SE_integ}")
        
    elif mode == "full":
        print('\nFULL ACF\n')
        autocorr_geq0=[x for x in autocorr if x >= 0]
        positive_integr_range=range(len(autocorr_geq0))
        autocorr_leq0=[-x for x in autocorr if x < 0]
        negative_integr_range=range(len(autocorr_geq0),len(time))
        t_auto_integ = integral.simps(autocorr_geq0, positive_integr_range) + integral.simps(autocorr_leq0, negative_integr_range)
        print(f'Autocorrelation time: {t_auto_integ} MC steps')

        N_ind_integ = t_sim/t_auto_integ
        SE_integ =  np.sqrt(np.var(entropy_mean)/N_ind_integ) # std deviation of the mean
        print(f"Number of independent samples: {N_ind_integ}\nStandard deviation of the mean: {SE_integ}")

    return t_auto_integ, N_ind_integ, SE_integ
    
#------------------------------------------------------------------------------------------------#

def block_analysis(array, chunk_length):
    """
    Script to calculate the BSE to measure correlation in an array.
    It takes the desired array as input and the chunk length.
    Returns the standard deviation array.
    """
    considered_length = array.shape[0]//chunk_length *chunk_length
    n_blocks = array.shape[0]//chunk_length
    
    block_array=array[:considered_length].copy().reshape(n_blocks, chunk_length)
    
    avg_block=np.average(block_array, axis=1)
    var_block_array=np.var(avg_block)/(avg_block.shape[0]-1)
    
    return np.sqrt(var_block_array)

#------------------------------------------------------------------------------------------------#

def block_analysis_estimation(block_results, time_sim):
    """
    Function taking as input the block analysis results and the total number of MC steps.
    Elaborates them to retrieve the autocorrelation time and the number of independent samples
    """
    n_ind = block_results[1]
    
    BSE_ind = block_results[0][n_ind-1]
    var_ind = np.var(block_results[0])

    corr_time = BSE_ind**2/var_ind * time_sim
    ind_samples = var_ind/BSE_ind**2
    print(f"autocorrelation time: {corr_time} MC steps\nN independent samples: {ind_samples}")

#------------------------------------------------------------------------------------------------#

def assessing_convergence(entropy_mean):
    """
    Function taking as input a pd Series (mean of entropy), converts it to numpy array 
    and performs block averaging on it.
    Returns the block entropy and the chunk length.
    """
    entropy_nparray=np.asarray(entropy_mean)
    entropy=entropy_mean.shape[0]//2
    block_entropy=np.zeros(entropy)
    chunk_length=2
    
    while chunk_length < entropy:
        block_entropy[chunk_length]=block_analysis(entropy_nparray, chunk_length)
        chunk_length+=1
        
    return block_entropy, chunk_length
    
#------------------------------------------------------------------------------------------------#
