import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import sys
from lib import statistics as EM_stat

data = np.random.normal(0, 1, 10000) # generate random normal dist m=0 sigma=1 1000points
                                    # porque con estos parametros funciona el test de Kolmogorov
#print(data)
mean, var = stats.distributions.norm.fit(data) # get mean and var
print(mean, var)


x = np.linspace(-5, 5, 100)
#print(x)

fitted_data = stats.distributions.norm.pdf(x, mean, var)

plt.hist(data, density=False, bins=21, align='left', color='blue', edgecolor='black', linewidth=1)
plt.plot(x, fitted_data,'r-')
plt.show()

l_kurtosis, l_skew, stat_shapiro_wilk, p_shapiro_wilk,\
stat_agostino_pearson, p_agostino_pearson, \
stat_kolmogorov_smirnov, p_kolmogorov_smirnov = EM_stat.does_data_follow_a_normal_distr(data, 1)  # 1: show results
