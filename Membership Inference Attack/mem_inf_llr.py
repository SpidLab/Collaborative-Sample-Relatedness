import numpy as np 
from scipy import stats
import pandas as pd
import math
import random


def compute_thresh(mix_stats):
    #input data
    population_gp_file = pd.read_csv("reference_pop_chrm2.csv", sep = '\t', index_col=0)
    # get population alternate allele frequencies
    pop_stats = (population_gp_file.apply(lambda x: (x==1).sum(), axis=1) + 2 * population_gp_file.apply(lambda x: (x==2).sum(), axis=1)) / (2.0 * 200)

    #control
    users = ['NA12286', 'HG01610', 'NA20815', 'HG00282', 'HG00375', 'NA12489', 'NA11995', 'NA20806', 'HG00267', 'HG00242', 'HG00176', 'HG00264', 'HG00102', 'HG01756', 'HG01522', 'HG00141', 'NA11843', 'HG01694', 'HG00371', 'HG00330', 'NA20521', 'HG00107', 'NA20520', 'NA12045', 'NA20524', 'HG00178', 'HG01771', 'NA12273', 'HG00237', 'HG01625', 'NA20768', 'NA20540', 'HG00119', 'HG01673', 'NA20811', 'HG00145', 'NA11829', 'HG00133', 'HG00123', 'HG01699']

    llr = []
    for i in range(len(users)):
        gt = population_gp_file[users[i]]

        comb_df = pd.concat([pop_stats, mix_stats, gt], axis = 1, join='inner', sort = False)
        comb_df.rename(columns={0:'pop_aaf'}, inplace=True)
        comb_df.rename(columns={'case_aaf':'mix_aaf'}, inplace=True)

        # compute llr
        comb_df['llr'] = comb_df.apply(lambda x: log_lr(x[str(users[i])], pd.to_numeric(x['mix_aaf']) , x['pop_aaf']), axis=1)
        ov_llr = comb_df['llr'].sum()
        llr.append(ov_llr)
        
    llr.sort(reverse = True)
    thresh = llr[2]
    return thresh

def log_lr(target, mixture, reference):
    target_norm = target/2.0
    llr = target_norm * math.log(mixture / reference) + (1 - target_norm) * math.log((1 - mixture) / (1 - reference))
    return llr


if __name__ == "__main__":

    # input data
    population_gp_file = pd.read_csv("reference_pop_chrm2.csv", sep = '\t', index_col=0)
    cols = population_gp_file.columns.values
    # get population alternate allele frequencies
    pop_stats = (population_gp_file.apply(lambda x: (x==1).sum(), axis=1) + 2 * population_gp_file.apply(lambda x: (x==2).sum(), axis=1)) / (2.0 * 200)
    # parameters
    alpha_value = 0.05 
    
    # get mixture aaf
    or_mixture_gp_file = pd.read_csv("mixture_chrm2.csv", sep = '\t', index_col=0)
    users = ['HG00129', 'HG00382', 'HG00244', 'NA20586', 'NA12249', 'HG00383', 'NA12760', 'HG00183', 'HG01525', 'NA20531', 'HG00174', 'NA12341', 'HG00325', 'HG00251', 'NA06994', 'HG00128', 'HG00235', 'HG01501', 'NA11933', 'HG00339', 'HG00281', 'NA20778', 'NA20804', 'HG01709', 'HG01770', 'NA20581', 'NA20512', 'NA12272', 'HG00266', 'HG00367', 'NA20762', 'HG00336', 'HG00099', 'HG01613', 'HG01779', 'HG00284', 'NA20805', 'NA12813', 'HG01791', 'HG00112']

    outf = open('llr_case_power_vs_num_snps.txt', 'a+')
    outf.write('num_snps,tp,power\n')
    for end_snp in range (10, 1001, 10):
        or_gwas = pd.read_csv("gwas_res_mixture.csv", sep = '\t', index_col=0)
        or_gwas = or_gwas.sort_values(by ='p_val', ascending=True)
        or_gwas = or_gwas[0:end_snp]

        # get mixture aaf
        mix_stats = or_gwas['case_aaf']
        thres = compute_thresh(mix_stats)
        print (thres)
        pred_y = []
        for i in range(len(users)):
            gt = or_mixture_gp_file[users[i]]

            comb_df = pd.concat([pop_stats, mix_stats, gt], axis = 1, join='inner', sort = False)
            comb_df.rename(columns={0:'pop_aaf'}, inplace=True)
            comb_df.rename(columns={'case_aaf':'mix_aaf'}, inplace=True)
            
            # compute llr
            comb_df['llr'] = comb_df.apply(lambda x: log_lr(x[str(users[i])], x['mix_aaf'] , x['pop_aaf']), axis=1)
            ov_llr = comb_df['llr'].sum()
            
            print(str(users[i]) + ',' + str(ov_llr))
            if ov_llr > thres:
                print ('member')
                pred_y.append(1)
            else:
                print ('non_member')
                pred_y.append(0)

        outf.write(str(len(comb_df)) + ',' + str(sum(pred_y)) + ',' + str(sum(pred_y) / len(users)) + '\n')

