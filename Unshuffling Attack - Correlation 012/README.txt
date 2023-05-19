1- compute_correlation.py
	-computes the correlation of SNPs between the same dataset
2- compute_correlation_2_datasets.py
	-computes the correlation of SNPs between two different datasets
	- changing the value of original_userIDs = userIDs[0:900], yields different results (900 out of 942 possible users)
3- unshuffling_attack.py
	- the first full implementation which checks the accuracy of unshuffling (same dataset, different number of users)
	- original_userIDs = userIDs[0:500] Modify these variables to get different results
	- attacker_userIDs = userIDs[0:500]
4- unshuffling_attack_noisy.py
	- checks unshuffling accuracy when we also add noise to one of the datasets
5- fast_unshuffling_attack_noisy.py
	- only computes the correlation tables for 1 SNP, therefore it is faster in running time
6- more_SNPs_unshuffling_attack.py
	- different from other scripts, the set of SNPs provided to the attackers ( in attacker_df ) is larger