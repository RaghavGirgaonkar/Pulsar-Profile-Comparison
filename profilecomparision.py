#!env python3

import numpy as np
import scipy.fft as fft
import lmfit
import statsmodels.stats.diagnostic as smdiag
import os
import pandas as pd
import seaborn as sns
# sns.set()

def fft_shift(prf, Δ): 
    """
        Circularly shift a profile prf by fractional rotation Δ.
        Δ ∈ [0,1]
        
        This uses the Fourier shift theorem.
        https://en.wikipedia.org/wiki/Fourier_transform#Translation_/_time_shifting
    """

    prf_fft = fft.rfft(prf) 
    fft_len = len(prf_fft) 
    prf_fft_shifted = np.exp(-1j*2*np.pi*np.arange(fft_len)*Δ) * prf_fft 
    prf_shifted = fft.irfft(prf_fft_shifted) 
    return prf_shifted 

def prf_shift_scale(prf, Δ, a, b): 
    """
        Circularly shift a profile prf by fractional rotation Δ, 
        scale it by amplitude a and add a baseline b.
    """    
    
    prf_shifted = fft_shift(prf, Δ) 
    return b + a*prf_shifted
    
def prf_diff_fn(prf1, prf2):
    def prf_diff(params):
        Δ = params['Δ']
        a = params['a']
        b = params['b']
        prf2_shifted_scaled = prf_shift_scale(prf2, Δ, a, b)
        return prf1 - prf2_shifted_scaled
        
    return prf_diff
    
def guess_align_params(prf1, prf2): 
    Δ0 = (np.argmax(prf1) - np.argmax(prf2))/len(prf2) 
    a0 = (np.max(prf1)-np.min(prf1)) / (np.max(prf2)-np.min(prf2))  
    b0 = np.mean(prf1)/a0 - np.mean(prf2)  
    return Δ0,a0,b0

def prf_align(prf1, prf2):
    """
        Align two profiles prf1 and prf2 using the Taylor method.
        Taylor 1992 (DOI: 10.1098/rsta.1992.0088)
    """

    Δ0,a0,b0 = guess_align_params(prf1, prf2)
    params = lmfit.Parameters()
    params.add('Δ', value=Δ0)
    params.add('a', value=a0)
    params.add('b', value=b0)
    
    prf_diff = prf_diff_fn(prf1, prf2)
    
    res = lmfit.minimize(prf_diff, params)
    
    Δ = res.params['Δ']
    a = res.params['a']
    b = res.params['b']
    
    return prf_shift_scale(prf2, Δ, a, b) 

def prf_compare(prf1, prf2):
    """
        Compare two profiles prf1 and prf2 using the Ljung-Box test.
        The outputs are the Ljung-Box statistic and the corresponding p-value.
        https://en.wikipedia.org/wiki/Ljung-Box_test
    """

    prf2a = prf_align(prf1, prf2)
    diff = prf1 - prf2a
    result = smdiag.acorr_ljungbox(diff, lags=len(diff)/2, return_df=True).iloc[0]
    return prf2a, diff, result

def get_chisqrs(prf,diff,nbins):
    """ Get the reduced chi sqrs using the RMS calculated by the off pulse region of the profile after scaling
        and aligning it to the template. Chisqrs = Sum(diff**2/op_rms)/(nbins -1)
        Off pulse region has been decided after manually examining template.
     """ 
    off_pulse = np.zeros(39)
    off_pulse[:20] = prf[:20]
    off_pulse[20:] = prf[45:] #Making off pulse region
    # print("Off pulse Region ",off_pulse)
    op_rms = np.var(off_pulse) #Rms
    # print("Off pulse RMS ",op_rms)
    s = 0
    for d in diff:
        s += d**2/op_rms

    s = s/(nbins - 1)
    # print("Chisqr value = ",s)

    return s 
    

if __name__=="__main__":

    import sys
    import matplotlib.pyplot as plt

    prf0file = sys.argv[1]
    prffiles = sys.argv[2:]
    
    prf0 = np.genfromtxt(prf0file)
    
    # f = open("Ljung-Box-Test-Results-New.txt",'a')
    # f.seek(0)
    # f.truncate()
    directory = os.getcwd()
    dirr = directory.split("/")
    f2 = open(dirr[-1]+"_Chisqrs.txt",'a')
    f2.seek(0)
    f2.truncate()
    Chi_sqrs = []
    difference_off_pulses = []
    all_differences = []
    prf_as = []
    for prf1file in prffiles:
        names = prf1file.split("_")
        print("Running prfcmp for", prf1file, ' and ',prf0file)
        # f.write("Running prfcmp for "+str(prf1file) + " and " + str(prf0file)+"\n")
        prf1 = np.genfromtxt(prf1file)
        
        if len(prf1) != len(prf0):
            print("[ERROR] Unequal Nbins found in {}. I will skip this file.\n".format(prf1file))
            continue
        
        prf1a, diff, result = prf_compare(prf0, prf1)
        print("Peak Difference = ",np.max(diff))
        # prf_as.append(prf1a)
        # all_differences.append(diff)

        # diff_off_pulse = np.zeros(39)
        # diff_off_pulse[:20] = diff[:20]
        # diff_off_pulse[20:] = diff[45:] #Making off pulse region
        # difference_off_pulses.append(diff_off_pulse)

        chisqr = get_chisqrs(prf1a,diff,np.shape(diff)[0])
        print("Chisqr value = ",chisqr)
        Chi_sqrs.append(chisqr)
        f2.write(names[1] + " " + str(chisqr)+"\n")

        f1 = open(names[1] +"-Differences.txt","w")
        # f1.seek(0)
        # f1.truncate()
        for d in diff:
            f1.write(str(d)+"\n")

        f1.close()
        
        # f.write("Ljung-Box test results\n----------------------\nLB statistic = "+ str(result['lb_stat'])+"\nLB p-value = "+str(result['lb_pvalue'])+"\n----------------------\n")
        # print("Ljung-Box test results")
        # print("----------------------")
        # print("LB statistic = ", result['lb_stat'])
        # print("LB p-value = ", result['lb_pvalue'])
        # print("----------------------\n")
        
        lin = np.linspace(0,1,64)
        # fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
        # fig  = plt.figure()
		
		
        plt.clf()
        # plt.subplot(211)
        ax1 = plt.subplot(211)
        # plt.title("Cycle 39 Band3")
        ax1.plot(lin,prf0, label="Reference Profile")
        ax1.plot(lin,prf1a, label="Scaled and Aligned Profile")
        plt.legend()
        # plt.subplot(212)
        ax2 = plt.subplot(212, sharex = ax1)
        ax2.plot(lin,diff, label="difference")
        ax2.set_xlabel("Phase")

        
        plt.legend(loc = 'best')
        plt.savefig(prf1file+".png")
        # plt.show()
    # np.savetxt("All-Differences.txt",all_differences)
    # np.savetxt("All-Offpulse-Diff.txt",difference_off_pulses)
    # prf_as_max = np.max(prf_as,axis = 0)
    # prf_as_min = np.min(prf_as,axis = 0)
    # np.savetxt("prf-as-maxmin.txt",np.transpose([prf_as_max,prf_as_min]))
    # ax3 = plt.subplot(313)
    # lin_new = np.linspace(0,len(Chi_sqrs),len(Chi_sqrs))
    # plt.figure()
    # plt.title("Chi Square Values by Epoch")
    # plt.scatter(lin_new,Chi_sqrs,label = "Chi sqrs for each epoch")
    # plt.axhline(y = 1,color= 'red')
    # plt.savefig("-Chisqrs.png")
    # f.close()
    f2.close()
    # lin = np.linspace(0,1,64)
    # plt.figure()
    # prf_res = []
    # for file in os.listdir(directory):
    #     print(file)
        
    #     if file.endswith("-Residuals.txt"):
    #         d = pd.read_csv(str(file),header = None)
    #         res = np.array(d[0])
    #         prf_res.append(res)
    #         plt.title(dirr[-1]+" Band 3 Differences")
    #         # plt.plot(lin,res)

    # prf_res = np.array(prf_res)
    # print(np.shape(prf_res))
    # # f3 = open("Differences-Stdvs.txt",'w')
    # res_medians = np.median(prf_res,axis = 0)
    # res_stds = np.std(prf_res,axis = 0)
    # results = np.transpose([res_medians,res_stds])
    # np.savetxt("Differences-STDVS.txt",results)
    # # plt.errorbar(lin,res_medians,res_stds)
    # plt.plot(lin,res_medians,color = 'blue')
    # plt.fill_between(lin,res_medians - res_stds,res_medians + res_stds,alpha = 0.5)
    # plt.xlabel("Phase")
    # # plt.savefig(dirr[-1] + " Subband 7 Differences.pdf") 
    # # plt.savefig(dirr[-1] + " Subband 7 Differences.pdf")       
    # plt.show()
        
    
