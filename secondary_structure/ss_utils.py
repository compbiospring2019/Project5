# Utils for Project 3
import os

# Get the parent directory of this code
this_script = os.path.abspath(__file__)
parent_directory = os.path.dirname(this_script)


def read_dist(file_path, dir=None):
    # Prior - probability of class label
    # Dists - list of means and standard deviations
    if dir:
        file_path = os.path.join(dir, file_path)
    with open(file_path, 'r') as f:
        prior = float(f.readline())
        dists = [[float(x) for x in line] for line in f]
    return prior, dists


def compute_percentages(ss_prediction):
    label_c_len = float(len([x for x in ss_prediction if x == 'C']))
    label_e_len = float(len([x for x in ss_prediction if x == 'E']))
    percents = {
        'C': label_c_len / len(ss_prediction),
        'E': label_e_len / len(ss_prediction)
    }
    percents['H'] = 1.0 - percents['C'] - percents['E']
    return percents
