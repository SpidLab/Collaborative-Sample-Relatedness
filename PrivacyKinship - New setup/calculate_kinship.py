from email import parser
import numpy as np 
import pandas as pd
import itertools
import math
import random
import editdistance
from itertools import permutations 
import argparse


def calculate_phi(df, user1, user2, userList):
    n11 = 0
    n02 = 0
    n20 = 0
    n_1 = 0
    n1_ = 0
    index1 = userList.index(user1)
    index2 = userList.index(user2)
    for row in df.itertuples(index=False):
        if row[index1] == 1 and row[index2] == 1:
            n11 += 1
        if row[index1] == 0 and row[index2] == 2:
            n02 += 1
        if row[index1] == 2 and row[index2] == 0:
            n20 += 1
        if row[index1] == 1:
            n1_ += 1
        if row[index2] == 1:
            n_1 += 1
    if n1_ != 0:
        phi = (2*n11-4*(n02+n20)-n_1+n1_)/(4*n1_)
    else:
        phi = -999
    return phi

def compute_coefficients_dictionary(url, nSNPs):
    threshold = 0.08
    TP = TN = FP = FN = 0
    df = pd.read_csv(url, sep=',', index_col=0)
    df = df.sample(nSNPs, random_state=88)
    userList = df.columns.to_list()
    related_df = pd.read_csv("Data/related_ids.csv", sep=',', index_col=0)
    unrelated_df = pd.read_csv("Data/unrelated_ids.csv", sep=',', index_col=0)
    related_coeff = {}
    unrelated_coeff = {}
    # the general dictionary that has the following stucture
    # {('10', '100000'): 0.25, ...}
    for i in range(30):
        first_user = related_df.at[i, "user1"]
        second_user = related_df.at[i, "user2"]
        phi_val_left = calculate_phi(df, str(first_user), str(second_user), userList)
        phi_val_right = calculate_phi(df, str(second_user), str(first_user), userList)
        phi_val = max(phi_val_left, phi_val_right) # keeping the largest \phi
        related_coeff[(first_user, second_user)] = phi_val

        first_user = unrelated_df.at[i, "user1"]
        second_user = unrelated_df.at[i, "user2"]
        phi_val_left = calculate_phi(df, str(first_user), str(second_user), userList)
        phi_val_right = calculate_phi(df, str(second_user), str(first_user), userList)
        phi_val = max(phi_val_left, phi_val_right)  # keeping the largest \phi
        unrelated_coeff[(first_user, second_user)] = phi_val

    for coeff in related_coeff.values():
        if coeff > threshold:
            TP += 1
        else:
            FN += 1
    for coeff in unrelated_coeff.values():
        if coeff > threshold:
            FP += 1
        else:
            TN += 1

    return TP, TN, FP, FN

if __name__ == "__main__":
    arr = [50, 250, 500, 1000, 2500]
    epsilon = [3, 4, 5]
    for eps in epsilon:
        print("Epsilon = " + str(eps))
        for nSNPs in arr:
            variables = compute_coefficients_dictionary(
                "Data/Eye_color_noisy/run_0/noisy_data_mixture_eps_" + str(eps) + ".0.csv", nSNPs)
            TP = variables[0]
            TN = variables[1]
            FP = variables[2]
            FN = variables[3]
            accuracy = (TP + TN) / (TP + TN + FP + FN)
            precision = TP / (TP + FP)
            recall = TP / (TP + FN)
            print("Number of SNPs: " + str(nSNPs) + "\tAccuracy: " + str(accuracy) + "\tRecall: " + str(
                recall) + "\tPrecision: " + str(precision))

    print("Epsilon = Inf")
    for nSNPs in arr:
        variables = compute_coefficients_dictionary("Data/Eye_color_relatives/mixture_with_all_relatives_3000_SNPs.csv",
                                                    nSNPs)
        TP = variables[0]
        TN = variables[1]
        FP = variables[2]
        FN = variables[3]
        accuracy = (TP + TN) / (TP + TN + FP + FN)
        precision = TP / (TP + FP)
        recall = TP / (TP + FN)
        print("Number of SNPs: " + str(nSNPs) + "\tAccuracy: " + str(accuracy) + "\tRecall: " + str(
            recall) + "\tPrecision: " + str(precision))

