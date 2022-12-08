# python3
# ################################################################## #
# files.py (C) March-2021 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# ------------------------------------------------------------------------
# Project: geneLength
#
# Purpose: module to deal with files
#
# ################################################################## #
import sys
sys.path.append('./lib/')

import pandas as pd
import collections

from scipy.stats import entropy
from scipy.stats import chi2_contingency, pearsonr, shapiro, normaltest, kstest, kurtosis, skew
import matplotlib.pyplot as plt


def estimate_shannon_entropy(dna_sequence):
    """
        Found in
        https://onestopdataanalysis.com/shannon-entropy/
        Calculate the Shannon entropy from a sequence

        Parameters
        ----------
        dna_sequence: str (*.tsv),
            The file (with path) where the data was presumably saved
        Output:
        ------
        entropy_value: float
            the Shannon entropy of the dna_sequence
     """

    bases = collections.Counter([tmp_base for tmp_base in dna_sequence])
    print(bases)
    print(bases.values())
    # define distribution
    dist = [x/sum(bases.values()) for x in bases.values()]

    # use scipy to calculate entropy
    return entropy(dist, base=2)


def pdf_of_normal_distr(x, mean, std):
    '''
        Calculate the "probability distribution function (pdf)" of the variable x
        given the next inputs:
            - x: input vector
            - mean: the mean of the normal distribution
            - std: the standard distribution of the normal distribution
            -output: also a vector with the pdf(x)
    '''

    pdf = np.exp(-(x-mean)**2/(2 * std**2))/(std * np.sqrt(2 * np.pi)) #prob. distr. function
    if 0:
            print("x:", x)
            print("pdf:", pdf)

    if 0: # calculate the integral if all the points (x) are at the same distance
        integral = 0
        step=x[1]-x[0]
        for p in pdf:
            integral += p*step
        print("\tpdf_of_normal_distr\n\t\tthe integral of pdf is:", integral)

    if 0: #plot the distribution
        plt.plot(x, step*pdf, linewidth=2, color='r')
        plt.show()
    return pdf

def does_data_follow_a_normal_distr(data, bool_show_results=1):
    '''
        Calculate the "probability distribution function (pdf)" of the variable x
        given the next inputs:
            - data: list
                input data
            - bool_show_results int: 0 or 1
                show the results
    '''

    l_kurtosis=kurtosis(data)
    l_skew=skew(data)
    # Shapiro-Wilk Test
    stat_shapiro_wilk,p_shapiro_wilk=shapiro(data)
    # D'Agostino and Pearson's Test
    stat_agostino_pearson, p_agostino_pearson = normaltest(data)
    #Kolmogorov Smirnov test
    stat_kolmogorov_smirnov, p_kolmogorov_smirnov = kstest(data, 'norm')

    if bool_show_results:
        ALPHA=0.05  # p_value threshold
        #
        # Kurtosis: valores atípicos para la normal (0: normal)
        print("Kurtosis:", l_kurtosis)
        # Skewness: hacia donde va la carga en la cola (0: no hay cola)
        print("Skewness:", l_skew)
        #Kolmogorov Smirnov test
        stat_kolmogorov_smirnov, p_kolmogorov_smirnov = kstest(data, 'norm')

        # Shapiro-Wilk Test
        print("Shapiro-Wilk Test")
        print('\tStatistics=%.3f, p_value=%.3f' % (stat_shapiro_wilk, p_shapiro_wilk))
        # For the sake of the interpretation
        if p_shapiro_wilk > ALPHA:
            print('\tSample looks Gaussian (fail to reject H0)')
        else:
            print('\tSample does not look Gaussian (reject H0)')

        # D'Agostino and Pearson's Test
        print("D'Agostino and Pearson's Test")
        print('\tStatistics=%.3f, p_value=%.3f' % (stat_agostino_pearson, p_agostino_pearson))
        # interpret
        if p_agostino_pearson > ALPHA:
            print('\tSample looks Gaussian (fail to reject H0)')
        else:
            print('\tSample does not look Gaussian (reject H0)')

        #Kolmogorov Smirnov test
        print("Kolmogorov Smirnov Test")
        print("\tks_statistic=", stat_kolmogorov_smirnov, ", p_value=", p_kolmogorov_smirnov)
         # interpret
        if p_kolmogorov_smirnov > ALPHA:
            print('\tSample looks Gaussian (fail to reject H0)')
        else:
            print('\tSample does not look Gaussian (reject H0)')

    return  l_kurtosis, l_skew, stat_shapiro_wilk, p_shapiro_wilk, stat_agostino_pearson, p_agostino_pearson, stat_kolmogorov_smirnov, p_kolmogorov_smirnov

'''
#test
df_tshirts = pd.DataFrame([ [48,63,33,47],
                            [35,246,42,27] ],
                          index=["Men","Women"],
                          columns=["Black","White","Red","Blue"])
print_results_chi2_contingency(df_tshirts, 1)
print_results_2var_correlation([48,63,33,47],  [35,246,42,27])
'''

if 0:
    seq = "ATGCGTCTCAAA"
    seq = "AAAATAAAAATAAAACG"
    print(entropy.__doc__)
    print (estimate_shannon_entropy(seq))