# functions used to process training data and make predictions about testing data
import secondary_structure.ss_utils as utils
from math import sqrt, exp
import os

empty_row = {'A': -1, 'C': -1, 'E': -1, 'D': -1, 'G': -1, 'I': -1, 'H': -1, 'K': -1, 'F': -1, 'M': -1, 'L': -1, 'N': -1, 'Q': -1, 'P': -1, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': -1, 'Y': -1}
acids_list = ['A', 'C', 'E', 'D', 'G', 'I', 'H', 'K', 'F', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y']

# Get the parent directory of this code
this_script = os.path.abspath(__file__)
parent_directory = os.path.dirname(this_script)


dists = {}


# max_prob - maximum probability the given feature values were observed given the specified class label
def maximum_likelihood(feature_values, dist_file, dir=parent_directory):
    dist_file = os.path.join(dir, dist_file)
    if dist_file not in dists:
        print('Reading in {}'.format(dist_file))
        dists[dist_file] = {'sigma': {}, 'mu': {}}
        with open(dist_file, 'r') as f:
            dists[dist_file]['prior'] = float(f.readline())
            for line_num in range(100):
                mean, std_dev = [float(x) for x in f.readline().split()]
                dists[dist_file]['sigma'][line_num] = std_dev
                dists[dist_file]['mu'][line_num] = mean

    prob = dists[dist_file]['prior']
    for feat_num in range(100):
        prob *= gnb(feature_values[feat_num], dists[dist_file]['mu'][feat_num], dists[dist_file]['sigma'][feat_num])
    return prob


def gnb(value, mean, std_dev):
    return exp(-0.5 * (value - mean) ** 2 / (std_dev ** 2)) / sqrt(2 * 3.14159 * (std_dev ** 2))


def classify(pssm_classify):
    predictions = []
    for row_num in range(len(pssm_classify)):
        # find feature values
        feature_values = []
        for row_offset in range(-2, 3):
            if row_num + row_offset < 0 or row_num + row_offset >= len(pssm_classify):
                # out of bounds
                feature_values.extend([-1] * 20)
            else:
                # not out of bounds
                row = pssm_classify[row_num + row_offset]
                feature_values.extend([row[k] for k in acids_list])
        # all feature values recorded
        # now find the maximum probability these features were observed given C, E, and H
        gnb_c = maximum_likelihood(feature_values, "C.dist")
        gnb_e = maximum_likelihood(feature_values, "E.dist")
        gnb_h = maximum_likelihood(feature_values, "H.dist")
        # prediction
        if max([gnb_c, gnb_e, gnb_h]) == gnb_c:
            prediction = 'C'
        elif max([gnb_c, gnb_e, gnb_h]) == gnb_e:
            prediction = 'E'
        else:
            prediction = 'H'
        predictions.append(prediction)

    # Compute percentages of each class label
    return utils.compute_percentages(''.join(predictions))
