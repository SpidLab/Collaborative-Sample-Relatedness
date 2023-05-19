import pandas as pd
import numpy as np
import scipy.stats as sc
from operator import itemgetter
import math
import json
import random

def compute_probabilities(first_list, second_list):
    dict = {
        '00': 0,
        '01': 0,
        '02': 0,
        '10': 0,
        '11': 0,
        '12': 0,
        '20': 0,
        '21': 0,
        '22': 0
    }
    prob = 1/len(first_list) #probability being added for each occurance of 00, 01 ....
    for i in range(len(first_list)):
        combination = str(first_list[i]) + str(second_list[i])
        dict[combination] += prob
    return dict

def compute_difference(first_dict, second_dict):
    diff = 0
    combination_list = ['00','01','02','10','11','12','20','21','22']
    for couple in combination_list:
        diff += abs(first_dict[couple] - second_dict[couple])
    return diff

# Read the dataframe that contains all the data
original_df = pd.read_csv("Data/Eye_color/original_3000_SNPs_mixture_chrm15.csv", sep=',', index_col=0)
attacker_df = pd.read_csv("Data/Eye_color/Noisy_data/run_0/noisy_data_mixture_eps_5.0.csv", sep=',', index_col=0)

#Select only the users we want to be part of computations
caseIDs = ['10', '10012', '10074', '1010', '10142', '1026', '1033', '1034', '1036', '1071', '11', '1114', '1121', '1125', '1142', '1163', '1173', '118', '1209', '1221', '1233', '1235', '126', '1281', '1288', '13', '1375', '1382', '1384', '1392', '14', '141', '1420', '144', '146', '1466', '1467', '1470', '1532', '1537', '154', '1566', '1627', '1668', '1678', '175', '1759', '1785', '1825', '187', '1877', '1879', '1885', '189', '1952', '1963', '1975', '1987', '1991', '1996', '2008', '2018', '2035', '2075', '2077', '2133', '2153', '2162', '2217', '2236', '2239', '2276', '2279', '2337', '2371', '2385', '2402', '2421', '2455', '2512', '2513', '2535', '2557', '2576', '258', '2616', '262', '2628', '2645', '266', '2675', '268', '2702', '2725', '2759', '276', '279', '2841', '2849', '285', '2858', '2887', '2901', '2907', '2949', '2960', '2964', '2966', '298', '2998', '3046', '3069', '311', '3130', '3140', '3196', '322', '323', '3257', '328', '33', '3335', '3375', '3381', '3394', '3397', '3399', '3468', '3486', '3499', '3547', '357', '3581', '3605', '3615', '3622', '3625', '3639', '366', '3689', '3695', '3766', '3780', '380', '3835', '3853', '3883', '3904', '3905', '3913', '3921', '3922', '3931', '3934', '3964', '3965', '397', '3973', '3979', '3992', '4006', '4018', '4029', '4032', '4034', '4048', '4066', '4084', '4127', '4171', '4185', '4230', '4346', '4381', '44', '440', '4411', '4428', '4490', '4507', '4569', '4577', '4580', '4605', '4632', '4643', '4646', '466', '4671', '4672', '4680', '4715', '4758', '481', '4813', '4814', '483', '4833', '4837', '4851', '4859', '488', '4897', '4943', '4990', '4998', '5105', '5146', '5167', '5182', '5185', '5195', '5207', '5263', '5288', '5302', '5310', '5332', '539', '5408', '5442', '5460', '549', '5515', '554', '5557', '556', '5581', '561', '563', '5675', '5719', '5723', '5745', '5803', '5830', '5833', '5881', '589', '5924', '5973', '5989', '6020', '6064', '6065', '6074', '6099', '610', '6104', '6120', '6124', '613', '6192', '62', '6249', '6356', '644', '6484', '6573', '6583', '6589', '6721', '6733', '6734', '675', '678', '6815', '6858', '6859', '6869', '6873', '688', '6882', '6936', '6981', '6993', '7020', '7047', '7055', '7080', '7092', '7144', '7153', '7156', '7195', '722', '726', '7285', '7312', '7358', '7360', '7379', '7405', '7413', '7446', '745', '7471', '7475', '748', '7485', '75', '7504', '754', '7547', '7549', '7565', '7585', '7610', '7611', '762', '763', '7759', '777', '7783', '7905', '792', '7929', '7950', '7954', '7955', '798', '7989', '799', '8', '803', '8046', '807', '8089', '813', '8134', '8138', '8155', '816', '8178', '8184', '820', '8203', '8266', '8291', '8314', '832', '8336', '8387', '8415', '8417', '8423', '845', '8522', '854', '8563', '8598', '862', '8635', '8737', '8787', '8788', '8803', '8808', '8818', '8841', '8852', '8874', '8890', '8988', '9078', '9133', '916', '9166', '918', '9185', '920', '922', '9220', '924', '9246', '926', '9262', '9375', '9399', '943', '9450', '9457', '951', '9516', '9527', '954', '9565', '9635', '966', '9697', '9702', '9720', '9726', '9743', '9747', '9799', '9803', '9812', '9839', '9869', '9969']
controlIDs = ['1', '10080', '10123', '1013', '1020', '1022', '1028', '1029', '1039', '1040', '1059', '1066', '1075', '1077', '1085', '1089', '1103', '1104', '1129', '1131', '1133', '1134', '1139', '1147', '1165', '1177', '1190', '1196', '123', '124', '1269', '1275', '1283', '1312', '1325', '1379', '1424', '1426', '1445', '1459', '1494', '1499', '1503', '1514', '1531', '1588', '159', '16', '160', '1600', '161', '168', '17', '1717', '1764', '177', '180', '1802', '1833', '1875', '1964', '1968', '1970', '1980', '200', '2000', '2004', '202', '2023', '2024', '2030', '204', '2056', '207', '2074', '2080', '2099', '210', '2116', '215', '2158', '216', '2166', '2177', '2183', '22', '2202', '2212', '2238', '2249', '2274', '2278', '2288', '2291', '2293', '2297', '2307', '2314', '2334', '2350', '2356', '2361', '2362', '2367', '241', '2412', '2453', '2454', '2456', '2457', '2498', '2507', '251', '252', '253', '2533', '2542', '2554', '2568', '2570', '26', '2633', '2648', '2660', '2663', '2681', '2687', '2700', '2707', '2718', '2720', '2739', '2764', '2797', '2826', '287', '2874', '288', '2881', '2886', '2890', '2920', '294', '295', '2958', '296', '2969', '2980', '2996', '3004', '3018', '3034', '305', '3078', '3090', '3104', '3167', '3186', '319', '325', '3265', '3275', '3280', '3289', '330', '3321', '3324', '3326', '337', '341', '3417', '3422', '3449', '345', '3454', '347', '3478', '349', '3494', '35', '352', '3531', '3532', '3543', '3569', '3582', '3598', '36', '3611', '363', '3642', '3643', '3644', '3655', '368', '3682', '3706', '3782', '38', '3821', '3834', '3847', '3865', '3881', '3882', '3908', '3930', '394', '3942', '3954', '3967', '40', '4009', '4030', '4036', '4038', '4055', '4056', '4057', '4080', '4095', '4106', '411', '4143', '4147', '4159', '4173', '4181', '4192', '42', '420', '4200', '429', '4349', '4356', '4359', '437', '439', '4399', '441', '4438', '4441', '4463', '4482', '4499', '4539', '4547', '4549', '4587', '4594', '4621', '463', '4648', '4679', '468', '4694', '4695', '4712', '473', '4756', '4797', '4801', '4825', '4879', '4932', '4935', '4941', '495', '4962', '497', '4983', '4985', '500', '502', '503', '5048', '5049', '5073', '5098', '5107', '5155', '5164', '5168', '5176', '5180', '5213', '5223', '5233', '5235', '5249', '5257', '5274', '5326', '533', '5342', '5397', '5425', '5432', '5447', '5456', '5475', '5478', '5487', '5495', '55', '5529', '5540', '5568', '5582', '5612', '5613', '5629', '567', '5683', '5701', '572', '574', '5741', '5743', '58', '5800', '581', '5820', '5821', '583', '5850', '5853', '5868', '5884', '5900', '5917', '595', '60', '6005', '6013', '602', '6035', '6045', '607', '609', '6103', '6115', '6187', '6191', '6206', '6238', '6239', '6244', '6248', '6269', '6284', '6285', '63', '6335', '6351', '636', '6363', '6367', '637', '6381', '64', '645', '646', '6483', '6491', '650', '651', '6545', '6560', '6617', '667', '669', '6690', '6699', '6706', '6731', '6758', '682', '684', '6851', '6861', '6866', '6889', '6910', '6928', '693', '6943', '6965', '6971', '6985', '6989', '7014', '704', '7053', '7070', '7072', '7082', '7089', '7098', '7114', '7155', '7163', '717', '7198', '7199', '721', '7214', '7218', '7232', '7252', '7266', '7280', '7292', '7320', '7322', '7393', '74', '7403', '7450', '7466', '749', '7510', '7520', '7568', '758', '761', '7636', '7639', '7657', '7685', '77', '7701', '775', '776', '7761', '7785', '7787', '779', '781', '782', '7841', '7876', '7915', '7923', '7945', '7967', '7981', '7992', '8058', '806', '8060', '8073', '8076', '808', '8090', '81', '8100', '811', '8129', '8135', '814', '8202', '8206', '8237', '8246', '8301', '8318', '836', '8363', '839', '8403', '8416', '842', '8449', '8463', '8476', '850', '8554', '8561', '865', '8675', '8701', '8731', '8796', '8864', '887', '8905', '8911', '8915', '8934', '8948', '8951', '8977', '898', '90', '903', '9072', '909', '9119', '912', '9126', '9135', '915', '9158', '9172', '9215', '9232', '925', '9276', '9323', '9324', '9351', '9367', '9372', '9380', '9421', '9427', '9430', '9486', '9489', '949', '9518', '9519', '952', '9541', '9561', '9568', '9589', '9616', '9672', '968', '972', '9749', '9781', '9796', '9871', '9876', '99', '990', '9920', '9928']
userIDs = caseIDs + controlIDs #holds all user IDs - 942
random.shuffle(userIDs)
original_userIDs = userIDs[0:500]
attacker_userIDs = userIDs[0:610]

#Select only the numerical values
original_df = original_df[original_userIDs]
attacker_df = attacker_df[attacker_userIDs]

#Select only the SNPs with unique MAF
allSNPs = attacker_df.index.values.tolist()
MAF = attacker_df.apply(lambda x: x.sum(), axis=1) / (2 * len(allSNPs))
MAF = MAF.drop_duplicates()
MAF = MAF.sort_values()
uniqueSNPlist = []
for index, value in MAF.items():
    uniqueSNPlist.append(index)

n = 200 # select only the top X SNPs, because of performance issues
SNPlist = random.sample(uniqueSNPlist, n)
# SNPlist = uniqueSNPlist[100:100+n]
attacker_df = attacker_df.loc[SNPlist]
original_df = original_df.loc[SNPlist]


#Initially compute the MAF values of the dataset provided by researchers (attacker_df) and the MAF that the attacker has (assuming original_df)
original_MAF = original_df.apply(lambda x: x.sum(), axis=1)/(2*len(original_userIDs))
attacker_MAF = attacker_df.apply(lambda x: x.sum(), axis=1) / (2 * len(attacker_userIDs))
print(original_MAF)
print(attacker_MAF)
MAF_difference = {}
for SNP1 in SNPlist:
    for SNP2 in SNPlist:
        MAF_difference[(SNP1, SNP2)] = abs(attacker_MAF[SNP1] - original_MAF[SNP2])
sorted_MAF_diff = dict(sorted(MAF_difference.items(), key=lambda item: item[1]))

with open("MAF.txt", 'w') as f:
    for key, value in sorted_MAF_diff.items():
        f.write('%s:%s\n' % (key, value))

#Suppose that attacker matches first_index with the actual SNP, and tries to find the correlated SNPs
first_index = list(sorted_MAF_diff.keys())[0][0]

original_correlation={}
#the general dictionary that has the following stucture
# {('rs12905389', 'rs6599770'): {'00': 0.53, '01': 0.19, '02': 0.02, '10': 0.073, '11': 0.13, '12': 0.003, '20': 0.01, '21': 0, '22': 0.023}, ...}
for j in range(len(SNPlist)):
    if first_index != SNPlist[j]:
        second_index = SNPlist[j]
        first_list_series = pd.Series(original_df.loc[first_index])
        second_list_series = pd.Series(original_df.loc[second_index])
        first_list = pd.Series(original_df.loc[first_index]).to_numpy()
        second_list = pd.Series(original_df.loc[second_index]).to_numpy()
        original_correlation[(first_index, second_index)] = compute_probabilities(first_list, second_list)

attacker_correlation={} #similar structure with original_correlation
for j in range(len(SNPlist)):
    if first_index != SNPlist[j]:
        second_index = SNPlist[j]
        first_list_series = pd.Series(attacker_df.loc[first_index])
        second_list_series = pd.Series(attacker_df.loc[second_index])
        first_list = pd.Series(attacker_df.loc[first_index]).to_numpy()
        second_list = pd.Series(attacker_df.loc[second_index]).to_numpy()
        attacker_correlation[(first_index, second_index)] = compute_probabilities(first_list, second_list)

index = 0
TP = 0
count1 = 0
for SNP in SNPlist:
    second_index = SNP
    index += 1
    if first_index != second_index:
        actual_probability_dictionary = original_correlation[(first_index, second_index)]
        min = 9999 #will be used to find the minimum difference between pairs of SNPs 3x3 matrices
        min_index = "rs..."
        other_diff = 0
        for i in range(len(SNPlist)):
            if SNPlist[i] != first_index:#we don't consider distance between the same SNP
                third_index = SNPlist[i]
                possible_probability_dictionary = attacker_correlation[(first_index, third_index)]
                diff = compute_difference(actual_probability_dictionary, possible_probability_dictionary)
                if diff < min:
                    min = diff
                    min_index = third_index
                elif diff == min:
                    if MAF_difference[(first_index, third_index)] < MAF_difference[(first_index, min_index)]:
                        min = diff
                        min_index = third_index

                if third_index == second_index:
                    other_diff = diff

        if min_index == second_index:
            # print(str(first_index) + "," + str(second_index) + "," + str(min_index) + "," + str(min)+","+'1')
            TP += 1
        else:
            print("---------------------------------")
            print(str(first_index) + "," + str(second_index) + "," + str(min_index) + "," + str(min) + "," + '0')
            print(actual_probability_dictionary)
            print(attacker_correlation[(first_index, second_index)])
            print(str(other_diff) + " " + str(sorted_MAF_diff[(first_index, second_index)]))
            print(attacker_correlation[(first_index, min_index)])
            print(str(min) + " " + str(sorted_MAF_diff[(first_index, min_index)]))
            print(index)
            if min != other_diff:
                print("Different")
            else:
                print("Same")
            print(index)
            count1 += 1
            print("--------------------------------")

print(TP)
print(count1)
print("Accuracy = " + str((TP+1)/n)) #since we already match the first one ourselves