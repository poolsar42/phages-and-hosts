#!/usr/bin/python3
import numpy as np
import os
import itertools

'''
array_1 = np.array([
    [0.733, 0.684, 0.520, 0.055, 0.519, 0.045, 0.966, 0.748, 0.514, 0.995, 0.522, 0.404, 0.906, 0.796, 0.438, 0.105],
    [0.130, 0.128, 0.565, 0.464, 0.046, 0.218, 0.923, 0.616, 0.417, 0.346, 0.765, 0.457, 0.852, 0.903, 0.116, 0.841],
    [0.842, 0.872, 0.877, 0.145, 0.897, 0.224, 0.694, 0.338, 0.054, 0.871, 0.006, 0.267, 0.421, 0.279, 0.494, 0.638],
    [0.960, 0.011, 0.035, 0.266, 0.335, 0.894, 0.062, 0.277, 0.553, 0.821, 0.733, 0.087, 0.957, 0.050, 0.179, 0.540]
])

array_1short = np.array([
    [0.519, 0.045, 0.966, 0.748, 0.514, 0.995, 0.522, 0.404, 0.906, 0.796, 0.438, 0.105],
    [0.046, 0.218, 0.923, 0.616, 0.417, 0.346, 0.765, 0.457, 0.852, 0.903, 0.116, 0.841],
    [0.897, 0.224, 0.694, 0.338, 0.054, 0.871, 0.006, 0.267, 0.421, 0.279, 0.494, 0.638],
    [0.335, 0.894, 0.062, 0.277, 0.553, 0.821, 0.733, 0.087, 0.957, 0.050, 0.179, 0.540]
])

array_2 = np.array([
    [0.484, 0.292, 0.966, 0.896, 0.213, 0.504, 0.269, 0.744, 0.686, 0.826, 0.521, 0.199, 0.026, 0.564, 0.717, 0.644],
    [0.785, 0.300, 0.221, 0.284, 0.768, 0.662, 0.710, 0.216, 0.731, 0.926, 0.063, 0.084, 0.242, 0.334, 0.962, 0.507],
    [0.319, 0.228, 0.289, 0.315, 0.703, 0.239, 0.315, 0.783, 0.239, 0.890, 0.958, 0.052, 0.519, 0.003, 0.721, 0.252],
    [0.988, 0.413, 0.853, 0.028, 0.988, 0.586, 0.445, 0.360, 0.624, 0.111, 0.779, 0.083, 0.539, 0.325, 0.092, 0.832]
])
'''
#Add some kind of transformation for compositional data
#We want to calculate correlations between the matrices, but they have different sizes
#rolling method might be the best

def corr2_coeff(A, B):
    # Stolen from https://stackoverflow.com/a/30143754
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]
    
    
    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)

    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None],ssB[None]))

def get_diagonals(corr_matrix):
    diagonal_list = []
    for i in range(0, corr_matrix.shape[1]):
        diagonal_list.append([i, corr_matrix.diagonal(offset=i)])
    for i in range(-1, (corr_matrix.shape[0] * -1), -1):
        diagonal_list.append([i, corr_matrix.diagonal(offset=i)])
    return diagonal_list

def calculate_corr_between_meme_matrices(array1, array2):
    max_value = 0.
    if (len(array1) < len(array2)):
        array1, array2 = array2, array1
    for i in range(len(array1) - len(array2) + 1):
        corr_matrix = corr2_coeff(array1.T, array2.T)
        corr_diagonals = get_diagonals(corr_matrix)
    
        #remove all diagonals smaller than 6 (just to have a threshold)
        corr_diagonals = [x for x in corr_diagonals if len(x[1]) >= 6]

        #Get the index of the maximum of the mean of the remaining diagonals
        corr_diagonals_mean = [np.mean(x[1]) for x in corr_diagonals]
        max_value = np.max(corr_diagonals_mean)
    #returns max correlation value, offset, length of the max correlation diagonal, array 1 length, array 2 length
    try:
        max_index = corr_diagonals_mean.index(max_value)
        return (max_value, corr_diagonals[max_index][0], len(corr_diagonals[max_index][1]), array1.shape[1], array2.shape[1])
    except ValueError:
        return(['NA', 'NA', 'NA', array1.shape[1], array2.shape[1]])
    
    


'''
np.savetxt('corr_matrix_test.tsv', corr_matrix, delimiter='\t')
'''

# motif_data = 'results_total.50+sites.tsv'
# motif_folder = '/home/asier/meme_suite_stuff/host_phage_results_total/'
# motif_result = 'results_total.50+sites.corr_probs.tsv'

# #Load motif info
# meme_motif_list = []
# with open(motif_data) as in_handle:
#     in_handle.readline()
#     for line in in_handle:
#         splitLine = line.rstrip('\n').split('\t')
#         meme_motif_list.append(splitLine[0])
#         print(repr(splitLine[0]))
# print(len(meme_motif_list))

# motif_info_result = open(motif_result, 'w')
# motif_info_result.write('\t'.join(['id_1', 'id_2', 'max_corr', 'offset', 'diag_len', 'id_1_len', 'id_2_len']) + '\n')

# for i, comb in enumerate(itertools.combinations_with_replacement(meme_motif_list, 2)):
#     print(i, comb)
#     array_1 = np.loadtxt(motif_folder + comb[0] + '.probs.tsv', delimiter = '\t')
#     array_2 = np.loadtxt(motif_folder + comb[1] + '.probs.tsv', delimiter = '\t')

#     results = list(map(str, calculate_corr_between_meme_matrices(array_1, array_2)))
#     print(results)
#     motif_info_result.write('\t'.join([comb[0], comb[1]] + results) + '\n')

# motif_info_result.close()


