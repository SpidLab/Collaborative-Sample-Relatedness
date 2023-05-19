import numpy as np
import csv
import math
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from pathlib import Path

def randomized_response(val, p, q):
    rand_val = np.random.uniform(0, 1)
    new_val = val
    if rand_val > p:
        if val == 0:
            new_val = 1
        elif val == 2:
            new_val = 1
        elif val == 1:
            if rand_val > p+q:
                new_val = 0
            else:
                new_val = 2
    return new_val

#add noise
for itr in range(5):
    for epsilon in range(30, 60, 10):
        eps = epsilon / 10.0
        print ('------', eps)
        p = np.exp(eps)/(np.exp(eps)+2)
        q = 1/(np.exp(eps)+2)

        print (p, q)
        data = pd.read_csv("Data/Eye_color_relatives/mixture_with_all_relatives_3000_SNPs.csv", sep = ",", index_col=0)
        user_id = data.columns
        for j in range(len(user_id)):
            data[str(user_id[j])] = data.apply(lambda x: randomized_response(x[str(user_id[j])],  p, q), axis=1)


        #the following code creates a directory if it doesn't exit. Otherwise it throws error
        output_file = 'noisy_data_mixture_eps_' + str(eps) + '.csv'
        output_dir = Path('Data/Eye_color_noisy/run_' + str(itr))
        output_dir.mkdir(parents=True, exist_ok=True)
        data.to_csv(output_dir / output_file)  # can join path elements with / operator